"""Construct assembly: Architecture -> intended amino-acid sequence + feature spans.

Rules enforced:
  * User-specified sequences are preserved exactly (never codon-"corrected" here).
  * Substitutions/deletions are applied only at validated positions.
  * Precursor vs mature chain is represented explicitly.
  * Components are NEVER shortened to meet a size target.
  * If any protein/domain component is unresolved, assembly reports that the
    intended protein is INCOMPLETE and must stay provisional (no fabrication).
"""
from __future__ import annotations

from pydantic import BaseModel

from .genetic_code import AA_ALPHABET
from .schemas import Architecture, Component


class FeatureSpan(BaseModel):
    component_index: int
    kind: str
    name: str
    start_aa: int  # 1-based inclusive, in the assembled precursor
    end_aa: int
    provenance: str


class AssemblyResult(BaseModel):
    precursor_protein: str  # full assembled chain (includes signal peptide if any)
    mature_protein: str  # chain after signal-peptide removal (== precursor if none)
    spans: list[FeatureSpan]
    is_precursor: bool
    complete: bool  # False if any protein/domain component unresolved
    incomplete_reason: str = ""
    applied_substitutions: list[str] = []
    junction_notes: list[str] = []


def assemble(arch: Architecture) -> AssemblyResult:
    order = arch.fusion_order or list(range(len(arch.components)))
    spans: list[FeatureSpan] = []
    pieces: list[str] = []
    cursor = 0  # 0-based offset
    incomplete = False
    incomplete_bits: list[str] = []
    signal_span: FeatureSpan | None = None

    for comp_i in order:
        comp = arch.components[comp_i]
        if comp.kind in ("protein", "domain") and (not comp.resolved or not comp.sequence):
            incomplete = True
            incomplete_bits.append(comp.name)
            # record a zero-length placeholder span so the UI can show the gap
            spans.append(FeatureSpan(
                component_index=comp_i, kind=comp.kind, name=comp.name,
                start_aa=cursor + 1, end_aa=cursor, provenance=comp.provenance.value,
            ))
            continue
        seq = comp.sequence.upper()
        start = cursor + 1
        end = cursor + len(seq)
        pieces.append(seq)
        cursor = end
        span = FeatureSpan(
            component_index=comp_i, kind=comp.kind, name=comp.name,
            start_aa=start, end_aa=end, provenance=comp.provenance.value,
        )
        spans.append(span)
        if comp.kind == "signal_peptide":
            signal_span = span

    precursor = "".join(pieces)

    # ---- substitutions / deletions (validated positions, 1-based) --------
    applied: list[str] = []
    if arch.substitutions and precursor:
        chars = list(precursor)
        # apply deletions last so earlier indices stay valid: sort by position desc
        for sub in sorted(arch.substitutions, key=lambda s: s.position, reverse=True):
            pos = sub.position
            if pos < 1 or pos > len(chars):
                raise ValueError(f"substitution position {pos} out of range 1..{len(chars)}")
            current = chars[pos - 1]
            if sub.from_aa and current != sub.from_aa.upper():
                raise ValueError(
                    f"substitution expected {sub.from_aa!r} at {pos} but found {current!r}"
                )
            if sub.to_aa == "":  # deletion
                del chars[pos - 1]
                applied.append(f"del{current}{pos}")
            else:
                to = sub.to_aa.upper()
                if to not in AA_ALPHABET:
                    raise ValueError(f"invalid substitution target {to!r}")
                chars[pos - 1] = to
                applied.append(f"{current}{pos}{to}")
        precursor = "".join(chars)
        # NOTE: substitutions can shift downstream spans; recompute conservatively
        spans = _recompute_spans_after_edit(arch, order, precursor)

    # ---- mature chain ----------------------------------------------------
    is_precursor = arch.is_precursor or (signal_span is not None)
    if signal_span is not None and precursor:
        mature = precursor[signal_span.end_aa :]
    else:
        mature = precursor

    # ---- junction notes --------------------------------------------------
    junctions: list[str] = []
    resolved_spans = [s for s in spans if s.end_aa >= s.start_aa]
    for a, b in zip(resolved_spans, resolved_spans[1:]):
        junctions.append(
            f"{a.name}(…{a.end_aa}) | {b.name}({b.start_aa}…) junction at residue {a.end_aa}/{b.start_aa}"
        )

    reason = ""
    if incomplete:
        reason = ("Intended protein is INCOMPLETE: unresolved component(s): "
                  + ", ".join(incomplete_bits)
                  + ". Design must remain provisional; no residues were invented.")

    return AssemblyResult(
        precursor_protein=precursor,
        mature_protein=mature,
        spans=spans,
        is_precursor=is_precursor,
        complete=not incomplete,
        incomplete_reason=reason,
        applied_substitutions=applied,
        junction_notes=junctions,
    )


def _recompute_spans_after_edit(
    arch: Architecture, order: list[int], edited_protein: str
) -> list[FeatureSpan]:
    """After substitutions/deletions, best-effort re-derive spans by re-walking
    component lengths. Deletions inside a component shrink it by one; we treat
    spans as contiguous in fusion order over the edited sequence length."""
    spans: list[FeatureSpan] = []
    cursor = 0
    # We cannot perfectly attribute a mid-component deletion without tracking the
    # edit offsets per component; for traceability we recompute using the ORIGINAL
    # component lengths and clamp to the edited length.
    for comp_i in order:
        comp = arch.components[comp_i]
        if comp.kind in ("protein", "domain") and (not comp.resolved or not comp.sequence):
            spans.append(FeatureSpan(
                component_index=comp_i, kind=comp.kind, name=comp.name,
                start_aa=cursor + 1, end_aa=cursor, provenance=comp.provenance.value,
            ))
            continue
        length = len(comp.sequence)
        start = cursor + 1
        end = min(cursor + length, len(edited_protein))
        spans.append(FeatureSpan(
            component_index=comp_i, kind=comp.kind, name=comp.name,
            start_aa=start, end_aa=end, provenance=comp.provenance.value,
        ))
        cursor = end
    return spans
