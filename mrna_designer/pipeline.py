"""Pipeline orchestration: wire the separated modules into one design run.

interpret -> (references already used in interpret) -> literature ->
assemble -> transcript -> optimize -> validate -> status.

Kept deliberately thin and typed; each stage is a module that can be tested in
isolation. This module makes NO scientific decisions of its own beyond sequencing
the stages and setting the final status from validation + resolution state.
"""
from __future__ import annotations

import hashlib
from datetime import date

from . import literature
from . import transcript as T
from .assemble import assemble
from .interpret import interpret
from .optimize import optimize
from .schemas import (
    Architecture,
    ConstructRequest,
    Design,
    DesignStatus,
    OptimizationObjectives,
    Provenance,
)
from .validate import get_size_profile, validate_design


def _design_id(req: ConstructRequest, arch: Architecture) -> str:
    h = hashlib.sha256()
    h.update((req.text or "").encode())
    h.update((req.protein_identifier or "").encode())
    h.update(arch.assembled_protein().encode())
    return "mrna_" + h.hexdigest()[:10]


def _topics_for(arch: Architecture, req: ConstructRequest) -> list[str]:
    topics = {"mrna_general", "utr", "polyA", "codon_optimization", "kozak"}
    if any(c.kind == "linker" for c in arch.components):
        topics.add("linker")
    if len([c for c in arch.components if c.kind in ("protein", "domain")]) >= 2:
        topics.add("fusion")
    return sorted(topics)


def run_pipeline(
    req: ConstructRequest,
    *,
    objectives: OptimizationObjectives | None = None,
    allow_network: bool = False,
    size_profile_id: str | None = None,
    skip_optimization: bool = False,
) -> tuple[Design, "ClarificationRound"]:
    from .schemas import ClarificationRound  # local import to avoid cycle in annotations

    # 1) interpret (also resolves references)
    arch, clar = interpret(req, allow_network=allow_network)
    design_id = _design_id(req, arch)

    # collect retrieval states from resolved components (reference fetches)
    retrieval_states = []
    # (component-level retrieval states are summarized; interpret attaches notes)

    # 2) literature (bounded)
    topics = _topics_for(arch, req)
    evidence, lit_state = literature.retrieve(topics, allow_network=allow_network)
    retrieval_states.append(lit_state)
    # attach which decision each piece supports (light mapping)
    for e in evidence:
        if "linker" in (e.title + e.finding).lower():
            e.supports_decision = "choice of fusion linker"
        elif "kozak" in (e.title + e.finding).lower():
            e.supports_decision = "Kozak translation-initiation context"
        elif "poly" in (e.title + e.finding).lower() or "utr" in (e.title + e.finding).lower():
            e.supports_decision = "UTR / poly(A) selection (context-dependent)"
        else:
            e.supports_decision = "general mRNA design background"

    is_demo = any(c.provenance == Provenance.FIXTURE for c in arch.components)

    design = Design(
        design_id=design_id, request=req, architecture=arch,
        evidence=evidence, retrieval_states=retrieval_states,
        created=date.today(), is_demo_fixture=is_demo,
    )

    # 3) assemble
    asm = assemble(arch)
    intended = asm.precursor_protein  # the CDS must encode the full precursor

    # If incomplete/blocked: stop before inventing a CDS.
    if not asm.complete or not intended:
        design.status = DesignStatus.BLOCKED
        design.validation = validate_design(
            intended_protein=intended or "X", assembly=asm,
            transcript=None, objectives=objectives,
        ) if intended else None
        return design, clar

    # 4) optimize (deterministic) unless caller only wants architecture
    host = (req.context.expression_host or (objectives.host if objectives else None) or "h_sapiens")
    obj = objectives or OptimizationObjectives(host=host)
    if objectives is None:
        obj.host = host

    if skip_optimization:
        design.status = DesignStatus.PROVISIONAL
        return design, clar

    opt = optimize(intended, obj)
    design.optimization = opt
    selected = next(c for c in opt.candidates if c.candidate_id == opt.selected_id)

    # 5) transcript
    # honor user-supplied UTRs if present
    five = None
    three = None
    if "five_utr" in req.preferred_elements:
        five = T.user_element("5utr", req.preferred_elements["five_utr"])
    if "three_utr" in req.preferred_elements:
        three = T.user_element("3utr", req.preferred_elements["three_utr"])
    tr = T.build_transcript(intended, selected.cds_rna, five_utr=five, three_utr=three)
    design.transcript = tr

    # 6) validate
    size_profile = get_size_profile(size_profile_id)
    report = validate_design(
        intended_protein=intended, assembly=asm, transcript=tr,
        objectives=obj, size_profile=size_profile,
    )
    design.validation = report

    # 7) status
    #   VALIDATED  : all hard checks pass AND no unresolved inputs/assumptions
    #   PROVISIONAL: otherwise (errors present, assumptions recorded, or demo data)
    #   BLOCKED    : reserved for unresolved essential inputs (handled earlier)
    clean = report.passed and not arch.unresolved and not arch.assumptions and not is_demo
    design.status = DesignStatus.VALIDATED if clean else DesignStatus.PROVISIONAL
    return design, clar
