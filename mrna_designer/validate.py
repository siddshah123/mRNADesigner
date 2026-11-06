"""Sequence and size validation.

Hard failures (severity=error):
  * invalid sequence alphabet
  * frameshift / incomplete coding triplets
  * internal stop codons in the intended CDS
  * translation != specified protein
  * incorrect fusion junctions
  * invalid feature coordinates
  * violations of explicit hard constraints (forbidden motifs, homopolymer cap)

Size reporting is explicit and separates the SOFTWARE input limit from
biological/delivery guidance (which is configurable, source-backed, and labeled
uncertain). Computational checks are never presented as experimental proof.
"""
from __future__ import annotations

from .assemble import AssemblyResult
from .genetic_code import (
    AA_ALPHABET,
    STOP,
    max_homopolymer_run,
    normalize_dna,
    translate,
)
from .schemas import (
    OptimizationObjectives,
    SizeReport,
    Transcript,
    ValidationItem,
    ValidationReport,
)

# Software input limit (NOT a biological recommendation). Generous ceiling for
# in-memory processing; distinct from delivery/biology guidance.
SOFTWARE_INPUT_LIMIT_NT = 50_000


def _rna_alphabet_ok(seq: str) -> bool:
    return all(c in "ACGU" for c in seq.upper())


def validate_design(
    *,
    intended_protein: str,
    assembly: AssemblyResult | None,
    transcript: Transcript | None,
    objectives: OptimizationObjectives | None,
    size_profile: "SizeProfile | None" = None,
) -> ValidationReport:
    report = ValidationReport()
    items = report.items  # append to the model's own list (Pydantic copies on init)

    # ---- protein alphabet -----------------------------------------------
    if not intended_protein or not all(c in AA_ALPHABET for c in intended_protein.upper()):
        items.append(ValidationItem(
            check="protein_alphabet", severity="error",
            message="Intended protein contains invalid amino-acid letters or is empty.",
        ))
        # cannot do much more without a protein
        report.translated_protein = ""
        return report
    items.append(ValidationItem(
        check="protein_alphabet", severity="pass",
        message=f"Protein alphabet valid ({len(intended_protein)} aa).",
    ))

    # ---- assembly completeness ------------------------------------------
    if assembly is not None and not assembly.complete:
        items.append(ValidationItem(
            check="assembly_complete", severity="error",
            message=assembly.incomplete_reason or "Assembly incomplete (unresolved component).",
        ))

    if transcript is None:
        items.append(ValidationItem(
            check="transcript_present", severity="unevaluated",
            message="No transcript assembled (design provisional/blocked); CDS checks not evaluated.",
        ))
        return report

    cds = normalize_dna(transcript.optimized_cds)  # DNA letters internally

    # ---- CDS alphabet ----------------------------------------------------
    if not _rna_alphabet_ok(transcript.optimized_cds):
        items.append(ValidationItem(
            check="cds_alphabet", severity="error",
            message="CDS contains letters outside {A,C,G,U}.",
        ))
    else:
        items.append(ValidationItem(check="cds_alphabet", severity="pass",
                                     message="CDS alphabet valid (A/C/G/U)."))

    # ---- frame -----------------------------------------------------------
    if len(cds) % 3 != 0:
        items.append(ValidationItem(
            check="reading_frame", severity="error",
            message=f"CDS length {len(cds)} nt is not a multiple of 3 (frameshift).",
        ))
        # translation check cannot proceed reliably
    else:
        items.append(ValidationItem(check="reading_frame", severity="pass",
                                     message="CDS length is a multiple of 3."))

        # ---- translation + internal stops -------------------------------
        try:
            full = translate(cds, to_stop=False)
        except ValueError as e:
            items.append(ValidationItem(check="translation", severity="error",
                                        message=f"Translation failed: {e}"))
            full = ""
        if full:
            # last residue should be stop; none before it
            body = full[:-1] if full.endswith(STOP) else full
            ends_with_stop = full.endswith(STOP)
            internal_stops = body.count(STOP)
            report.translated_protein = body
            if internal_stops > 0:
                items.append(ValidationItem(
                    check="internal_stop", severity="error",
                    message=f"{internal_stops} internal stop codon(s) found in the CDS.",
                ))
            else:
                items.append(ValidationItem(check="internal_stop", severity="pass",
                                             message="No internal stop codons."))
            if not ends_with_stop:
                items.append(ValidationItem(
                    check="terminal_stop", severity="error",
                    message="CDS does not end in a stop codon.",
                ))
            else:
                items.append(ValidationItem(check="terminal_stop", severity="pass",
                                             message="CDS ends in a stop codon."))
            # translation identity vs intended protein
            if body != intended_protein.upper():
                items.append(ValidationItem(
                    check="translation_identity", severity="error",
                    message="Translated CDS does NOT match the intended protein "
                    f"(translated {len(body)} aa vs intended {len(intended_protein)} aa).",
                ))
            else:
                items.append(ValidationItem(
                    check="translation_identity", severity="pass",
                    message="Translated CDS matches the intended protein exactly.",
                ))

    # ---- fusion junctions (coordinates valid + contiguous) --------------
    if assembly is not None:
        resolved = [s for s in assembly.spans if s.end_aa >= s.start_aa]
        ok = True
        prev_end = 0
        for s in resolved:
            if s.start_aa != prev_end + 1 or s.end_aa < s.start_aa:
                ok = False
                break
            prev_end = s.end_aa
        if resolved and prev_end != len(assembly.precursor_protein):
            ok = False
        items.append(ValidationItem(
            check="fusion_junctions",
            severity="pass" if ok else "error",
            message="Component spans are contiguous and cover the assembled protein."
            if ok else "Component spans are non-contiguous or mis-coordinated.",
        ))

    # ---- explicit hard constraints --------------------------------------
    if objectives is not None:
        hp = max_homopolymer_run(cds)
        if hp > objectives.max_homopolymer:
            items.append(ValidationItem(
                check="homopolymer_constraint", severity="error",
                message=f"Max homopolymer run {hp} exceeds hard cap {objectives.max_homopolymer}.",
            ))
        else:
            items.append(ValidationItem(check="homopolymer_constraint", severity="pass",
                                         message=f"Max homopolymer run {hp} within cap."))
        # forbidden motifs in the CDS
        from .optimize import _count_motifs

        fm = _count_motifs(cds, objectives.forbidden_motifs)
        if fm > 0:
            items.append(ValidationItem(
                check="forbidden_motifs", severity="error",
                message=f"{fm} forbidden-motif occurrence(s) remain in the CDS.",
            ))
        elif objectives.forbidden_motifs:
            items.append(ValidationItem(check="forbidden_motifs", severity="pass",
                                         message="No forbidden motifs present."))

    # ---- RNA structure (unevaluated unless a tool is wired) -------------
    items.append(ValidationItem(
        check="rna_structure", severity="unevaluated",
        message="RNA secondary-structure metrics not evaluated (no folding tool "
        "configured in this build). Not a pass; shown as unevaluated.",
    ))

    # ---- size report -----------------------------------------------------
    report.size = _size_report(intended_protein, transcript, assembly, size_profile)

    # size vs software limit (hard) ; biological guidance (soft, labeled)
    if report.size.total_rna_len_nt > report.size.software_input_limit_nt:
        items.append(ValidationItem(
            check="software_size_limit", severity="error",
            message=f"Total RNA {report.size.total_rna_len_nt} nt exceeds the SOFTWARE "
            f"input limit {report.size.software_input_limit_nt} nt. Component contributions: "
            + ", ".join(f"{k}={v}" for k, v in report.size.component_contributions.items()),
        ))
    else:
        items.append(ValidationItem(check="software_size_limit", severity="pass",
                                     message=f"Total RNA {report.size.total_rna_len_nt} nt within software limit."))

    # biological size guidance (SOFT, labeled uncertain — never a hard failure)
    if size_profile is not None and size_profile.profile_id != "none":
        if report.size.total_rna_len_nt > size_profile.soft_max_nt:
            items.append(ValidationItem(
                check="biological_size_guidance", severity="warning",
                message=f"Total RNA {report.size.total_rna_len_nt} nt exceeds the "
                f"{size_profile.name} soft ceiling (~{size_profile.soft_max_nt} nt). "
                f"{size_profile.guidance} This guidance is UNCERTAIN "
                f"[source={size_profile.source_id}] and is NOT a hard failure. "
                "Component contributions: "
                + ", ".join(f"{k}={v}" for k, v in report.size.component_contributions.items()),
            ))
        else:
            items.append(ValidationItem(
                check="biological_size_guidance", severity="pass",
                message=f"Total RNA within the {size_profile.name} soft ceiling "
                f"(~{size_profile.soft_max_nt} nt); guidance is uncertain.",
            ))

    # computational-only disclaimer
    items.append(ValidationItem(
        check="scope_disclaimer", severity="unevaluated",
        message="These are COMPUTATIONAL checks only. They are not evidence of expression, "
        "secretion, activity, stability, or delivery.",
    ))

    return report


def _size_report(
    intended_protein: str,
    transcript: Transcript,
    assembly: AssemblyResult | None,
    size_profile: "SizeProfile | None",
) -> SizeReport:
    cds = normalize_dna(transcript.optimized_cds)
    five = normalize_dna(transcript.five_utr.sequence) if transcript.five_utr else ""
    koz = normalize_dna(transcript.kozak.sequence) if transcript.kozak else ""
    three = normalize_dna(transcript.three_utr.sequence) if transcript.three_utr else ""
    polyA_encoded = transcript.polyA.length if transcript.polyA.mode == "encoded" else 0
    total = len(transcript.mature_rna())

    contributions = {
        "5utr": len(five),
        "kozak": len(koz),
        "cds": len(cds),
        "3utr": len(three),
        "polyA_encoded": polyA_encoded,
    }

    guidance = "No biological size guidance profile selected; software limit only."
    guidance_src = None
    uncertain = True
    if size_profile is not None:
        guidance = size_profile.guidance
        guidance_src = size_profile.source_id
        uncertain = size_profile.uncertain

    return SizeReport(
        protein_len_aa=len(intended_protein),
        cds_len_nt=len(cds),
        cds_includes_stop=True,
        five_utr_len_nt=len(five),
        three_utr_len_nt=len(three),
        polyA_encoded_len_nt=polyA_encoded,
        polyA_mode=transcript.polyA.mode,
        total_rna_len_nt=total,
        software_input_limit_nt=SOFTWARE_INPUT_LIMIT_NT,
        biological_guidance=guidance,
        guidance_source_id=guidance_src,
        guidance_is_uncertain=uncertain,
        component_contributions=contributions,
    )


# --------------------------------------------------------------------------
# Configurable, source-backed biological size guidance (NOT a universal optimum)
# --------------------------------------------------------------------------
from pydantic import BaseModel  # noqa: E402


class SizeProfile(BaseModel):
    profile_id: str
    name: str
    soft_max_nt: int
    guidance: str
    source_id: str
    uncertain: bool = True


SIZE_PROFILES: dict[str, SizeProfile] = {
    "none": SizeProfile(
        profile_id="none", name="No biological guidance",
        soft_max_nt=10**9,
        guidance="No biological size constraint applied.",
        source_id="n/a", uncertain=True,
    ),
    "lnp_ivt_mrna": SizeProfile(
        profile_id="lnp_ivt_mrna", name="LNP-delivered IVT mRNA (illustrative)",
        soft_max_nt=5000,
        guidance="Illustrative soft ceiling (~5 kb) often cited for efficient LNP "
        "encapsulation/translation of IVT mRNA; context- and formulation-dependent.",
        source_id="lit_pardi_2018_review", uncertain=True,
    ),
}


def get_size_profile(profile_id: str | None) -> SizeProfile | None:
    if not profile_id:
        return None
    return SIZE_PROFILES.get(profile_id)
