"""Request interpretation: ConstructRequest -> provisional Architecture.

Deterministic, rule-based resolution of components (no fabrication). An LLM
could be slotted behind `interpret()` to parse free text more flexibly; the
contract is that it may only propose *structured* components/assumptions that
the deterministic pipeline then resolves against real sequences. The LLM never
emits nucleotides and never invents a protein sequence.

Resolution sources, in order of preference:
  1. User-pasted amino-acid sequence  (provenance=user, resolved)
  2. UniProt accession via live fetch  (provenance=reference, resolved if ok)
  3. Bundled synthetic demo protein    (provenance=fixture, resolved, flagged)
  4. Named real protein, not resolvable here -> UNRESOLVED dependency marked
"""
from __future__ import annotations

import re

from . import references as R
from .genetic_code import is_valid_protein, normalize_dna
from .schemas import (
    Architecture,
    Assumption,
    ClarificationOption,
    ClarificationQuestion,
    ClarificationRound,
    Component,
    ConstructRequest,
    Provenance,
)

# tag name aliases -> canonical library name
_TAG_ALIASES = {
    "his": "His6", "his6": "His6", "6xhis": "His6", "hexahistidine": "His6",
    "his-tag": "His6", "histag": "His6", "polyhistidine": "His6",
    "flag": "FLAG", "ha": "HA", "myc": "Myc", "c-myc": "Myc",
    "strep": "Strep-II", "strep-ii": "Strep-II", "strepii": "Strep-II",
    "v5": "V5",
}


def _detect_placement(text: str, keyword: str) -> str | None:
    """Look for 'N-terminal <keyword>' / 'C-terminal <keyword>' style cues."""
    t = text.lower()
    # within a small window around the keyword
    idx = t.find(keyword.lower())
    if idx == -1:
        return None
    window = t[max(0, idx - 30) : idx + 30]
    if re.search(r"\bn[- ]?term", window):
        return "N"
    if re.search(r"\bc[- ]?term", window):
        return "C"
    return None


def _resolve_protein_token(
    token: str, species: str | None, allow_network: bool
) -> list[Component]:
    """Resolve a single protein reference token into one or more Components.

    A precursor with a signal peptide is split into a `signal_peptide` component
    plus a mature `protein` component so the precursor/mature distinction is
    explicit downstream."""
    token = token.strip()
    # demo fixture?
    demo = R.get_demo_protein(token)
    if demo:
        ref = R.demo_as_reference(token)
        if ref.is_precursor and ref.signal_peptide and ref.mature_sequence:
            return [
                Component(
                    kind="signal_peptide", name=f"{ref.name} signal peptide",
                    sequence=ref.signal_peptide, species=ref.species,
                    accession=ref.accession, placement="N",
                    provenance=Provenance.FIXTURE, resolved=True,
                    source_id=ref.accession,
                    note="Signal peptide (precursor); " + ref.mature_note,
                ),
                Component(
                    kind="protein", name=f"{ref.name} mature chain",
                    sequence=ref.mature_sequence, species=ref.species,
                    accession=ref.accession, provenance=Provenance.FIXTURE,
                    resolved=True, source_id=ref.accession,
                    note="Mature chain of a BUNDLED SYNTHETIC DEMO FIXTURE.",
                ),
            ]
        return [Component(
            kind="protein", name=ref.name, sequence=ref.sequence,
            species=ref.species, accession=ref.accession,
            provenance=Provenance.FIXTURE, resolved=True,
            source_id=ref.accession,
            note="BUNDLED SYNTHETIC DEMO FIXTURE (not a real protein).",
        )]
    # uniprot accession?
    if R.looks_like_uniprot_accession(token):
        ref = R.fetch_uniprot(token, allow_network=allow_network)
        if ref.retrieval and ref.retrieval.ok and ref.sequence:
            return [Component(
                kind="protein", name=ref.name, sequence=ref.sequence,
                species=ref.species, accession=ref.accession,
                accession_version=ref.accession_version,
                provenance=Provenance.REFERENCE, resolved=True,
                source_id=ref.accession, note=ref.source,
            )]
        detail = ref.retrieval.detail if ref.retrieval else "no retrieval"
        return [Component(
            kind="protein", name=token, sequence="", accession=token,
            provenance=Provenance.ASSUMPTION, resolved=False, source_id=token,
            note=f"UniProt fetch did not resolve ({detail}). Sequence unresolved; "
            "paste the amino-acid sequence or enable network access.",
        )]
    # a plain protein NAME: try live UniProt name-search (research, no LLM)
    if allow_network:
        ref = R.search_uniprot(token, species=species, allow_network=True)
        if ref.retrieval and ref.retrieval.ok and ref.sequence:
            return [Component(
                kind="protein", name=ref.name, sequence=ref.sequence,
                species=ref.species, accession=ref.accession,
                accession_version=ref.accession_version,
                provenance=Provenance.REFERENCE, resolved=True,
                source_id=ref.accession, note=ref.source,
            )]
        detail = ref.retrieval.detail if ref.retrieval else "search failed"
        return [Component(
            kind="protein", name=token, sequence="", species=species,
            provenance=Provenance.ASSUMPTION, resolved=False,
            note=f"UniProt name-search did not resolve {token!r} ({detail}). "
            "Provide an accession or paste the sequence. No sequence was invented.",
        )]
    # offline: cannot resolve a bare name -> unresolved dependency (never invented)
    return [Component(
        kind="protein", name=token, sequence="", species=species,
        provenance=Provenance.ASSUMPTION, resolved=False,
        note="Named protein not resolved offline. Turn on 'Network retrieval' to "
        "look it up in UniProt automatically, give a UniProt accession, or paste "
        "the sequence. No sequence was invented.",
    )]


def _split_fusion(text: str) -> list[str]:
    """Very small heuristic splitter for 'A fused to B', 'A-B fusion', 'A + B'."""
    t = text
    for sep in [" fused to ", " fusion with ", " fused with ", " linked to ",
                " conjugated to ", " displayed on ", " and ", " + ", "/"]:
        if sep in t.lower():
            # split case-insensitively on first occurrence
            i = t.lower().find(sep)
            left, right = t[:i].strip(), t[i + len(sep):].strip()
            # strip a trailing/leading "fusion" noun from either side
            for junk in ("fusion", "construct", "a ", "an "):
                if right.lower().startswith(junk):
                    right = right[len(junk):].strip()
                if left.lower().endswith(" " + junk):
                    left = left[: -len(junk)].strip()
            return [left, right]
    return [text.strip()]


def interpret(req: ConstructRequest, *, allow_network: bool = False) -> tuple[Architecture, ClarificationRound]:
    """Produce a provisional Architecture and a (<=2 question) clarification round."""
    # If the user supplied an edited architecture directly, honor it verbatim.
    if req.architecture_override is not None:
        arch = req.architecture_override
        return arch, ClarificationRound(questions=[])

    text = req.text or ""
    species = req.context.protein_species
    components: list[Component] = []
    assumptions: list[Assumption] = []
    unresolved: list[str] = []

    # ---- core protein(s) -------------------------------------------------
    core_tokens: list[str] = []
    if req.pasted_protein_sequence:
        seq = re.sub(r"[^A-Za-z]", "", req.pasted_protein_sequence).upper()
        if is_valid_protein(seq):
            components.append(Component(
                kind="protein", name="User-supplied protein", sequence=seq,
                species=species, provenance=Provenance.USER, resolved=True,
                note="Amino-acid sequence supplied by user; preserved exactly.",
            ))
        else:
            unresolved.append("Pasted protein sequence contains invalid amino-acid letters.")
    elif req.protein_identifier:
        core_tokens = _split_fusion(req.protein_identifier)
    elif text:
        # extract a protein-ish phrase from the free text (very light heuristic)
        core_tokens = _split_fusion(text)

    for tok in core_tokens:
        for comp in _resolve_protein_token(tok, species, allow_network):
            components.append(comp)
            if comp.kind in ("protein", "domain") and not comp.resolved:
                unresolved.append(f"Unresolved protein identity/sequence: {comp.name!r}.")

    # ---- tags ------------------------------------------------------------
    blob = " ".join(filter(None, [text, req.protein_identifier or ""]))
    for alias, canonical in _TAG_ALIASES.items():
        if re.search(rf"\b{re.escape(alias)}\b", blob, re.IGNORECASE):
            tag = R.get_tag(canonical)
            if not tag:
                continue
            placement = _detect_placement(blob, alias)
            if placement is None:
                placement = "C"
                assumptions.append(Assumption(
                    field=f"tag:{canonical}:placement", value="C-terminal",
                    rationale="Terminus not specified; defaulted to C-terminal. Editable.",
                ))
            components.append(Component(
                kind="tag", name=canonical, sequence=tag["sequence"],
                placement=placement, provenance=Provenance.REFERENCE,
                resolved=True, source_id=f"tag:{canonical}", note=tag["source"],
            ))
            break  # one tag per pass (keep minimal; users can add more by editing)

    # ---- linker ----------------------------------------------------------
    protein_like = [c for c in components if c.kind in ("protein", "domain")]
    wants_linker = bool(re.search(r"\blinker\b", blob, re.IGNORECASE))
    if (len(protein_like) >= 2) and (wants_linker or " fused" in blob.lower()):
        lk = R.get_linker("GS_flexible_(GGGGS)x3")
        components.append(Component(
            kind="linker", name=lk["name"], sequence=lk["sequence"],
            placement="internal", provenance=Provenance.REFERENCE, resolved=True,
            source_id="linker:GS_flexible", note=lk["source"],
        ))
        assumptions.append(Assumption(
            field="linker", value=lk["name"],
            rationale="Two domains detected; defaulted to a flexible (GGGGS)x3 linker. Editable.",
        ))

    # ---- fusion order (N->C): N-tags, proteins/linkers in text order, C-tags
    order = _build_fusion_order(components)

    arch = Architecture(
        title=_title_from(req),
        components=components,
        fusion_order=order,
        substitutions=[],
        is_precursor=any(c.kind == "signal_peptide" for c in components),
        unresolved=unresolved,
        assumptions=assumptions,
        supported=True,
    )

    clar = _clarifications(req, arch)
    return arch, clar


def _title_from(req: ConstructRequest) -> str:
    if req.protein_identifier:
        return req.protein_identifier[:80]
    if req.text:
        return req.text[:80]
    return "Untitled construct"


def _build_fusion_order(components: list[Component]) -> list[int]:
    n_tags = [i for i, c in enumerate(components) if c.kind == "tag" and c.placement == "N"]
    c_tags = [i for i, c in enumerate(components) if c.kind == "tag" and c.placement == "C"]
    middle = [i for i, c in enumerate(components)
              if not (c.kind == "tag" and c.placement in ("N", "C"))]
    return n_tags + middle + c_tags


def _clarifications(req: ConstructRequest, arch: Architecture) -> ClarificationRound:
    """Prioritized, at most two questions, one round."""
    questions: list[ClarificationQuestion] = []

    # 1) unresolved protein identity (highest priority)
    unresolved_proteins = [c for c in arch.components
                           if c.kind in ("protein", "domain") and not c.resolved]
    if unresolved_proteins:
        c = unresolved_proteins[0]
        questions.append(ClarificationQuestion(
            key="protein_identity",
            question=f"I could not resolve {c.name!r} to an authoritative sequence. "
            "Provide a UniProt accession (e.g. P60568) or paste the amino-acid sequence?",
            options=[],
            impact="Without an authoritative sequence the coding region cannot be built; "
            "the design stays provisional.",
        ))

    # 2) species (encoded protein) vs expression host clarity
    if len(questions) < 2 and req.context.protein_species is None and (req.text or req.protein_identifier):
        questions.append(ClarificationQuestion(
            key="species_and_host",
            question="Which species is the ENCODED protein from, and which host/cell will "
            "EXPRESS it? (These are independent; a human protein need not be optimized for human cells.)",
            options=[
                ClarificationOption(label="Human protein, human cells", value="human/h_sapiens"),
                ClarificationOption(label="Human protein, mouse cells", value="human/m_musculus"),
                ClarificationOption(label="Human protein, E. coli", value="human/e_coli"),
            ],
            impact="Encoded species sets identity; host sets the codon-usage table used for optimization.",
        ))

    # 3) fusion orientation (if two protein domains and order ambiguous)
    protein_like = [c for c in arch.components if c.kind in ("protein", "domain")]
    if len(questions) < 2 and len(protein_like) >= 2:
        a, b = protein_like[0].name, protein_like[1].name
        questions.append(ClarificationQuestion(
            key="fusion_orientation",
            question=f"Fusion orientation (N->C): {a} then {b}, or {b} then {a}?",
            options=[
                ClarificationOption(label=f"{a}–{b}", value="as_written"),
                ClarificationOption(label=f"{b}–{a}", value="reversed"),
            ],
            impact="Orientation changes junctions, folding, and which terminus is free.",
        ))

    # 4) localization (only if still room and secretion is plausible/unspecified)
    if len(questions) < 2 and req.context.localization is None:
        questions.append(ClarificationQuestion(
            key="localization",
            question="Intended localization: secreted, membrane, or intracellular? "
            "(Determines whether a signal peptide is needed.)",
            options=[
                ClarificationOption(label="Secreted", value="secreted"),
                ClarificationOption(label="Intracellular", value="cytoplasmic"),
                ClarificationOption(label="Membrane", value="membrane"),
            ],
            impact="Secretion/membrane targeting requires a signal peptide, changing the precursor.",
        ))

    return ClarificationRound(questions=questions[:2])
