"""Deterministic, seeded, multi-objective codon optimization.

Strategy (documented, not a black box):
  1. Candidate generation
     - Baseline candidate: highest-usage codon at every position (max codon
       adaptation). This is the reference point, NOT assumed to be best.
     - Alternative candidates: seeded weighted sampling from the synonymous
       codon distribution (Random(seed+i)), giving diversity in GC/structure.
  2. Repair pass (applied to every candidate), only ever swapping a codon for
     a SYNONYMOUS one, so the encoded protein is invariant by construction:
     - break homopolymer runs longer than the limit,
     - remove user-forbidden motifs and restriction sites,
     - pull local-GC windows back toward the target band.
  3. Scoring: a weighted penalty (lower = better) over codon usage (1-CAI),
     global+local GC deviation, homopolymer excess, and motif/site hits.
  4. Selection: the minimum-penalty candidate; rationale records the tradeoff.

Reproducibility: identical (protein, objectives, seed) reproduce byte-identical
output. Version + seed + codon-table version are recorded in the result.
"""
from __future__ import annotations

import random

from .codon_tables import CodonTable, load_table
from .genetic_code import (
    SYNONYMS,
    gc_content,
    local_gc_extremes,
    max_homopolymer_run,
    normalize_dna,
    translate,
)
from .schemas import (
    OptimizationCandidate,
    OptimizationObjectives,
    OptimizationResult,
    SequenceMetrics,
)

OPTIMIZER_VERSION = "mrna-codon-opt/0.1.0"


def _count_motifs(dna: str, motifs: list[str]) -> int:
    s = normalize_dna(dna)
    total = 0
    for m in motifs:
        mm = normalize_dna(m)
        if not mm:
            continue
        start = 0
        while True:
            idx = s.find(mm, start)
            if idx == -1:
                break
            total += 1
            start = idx + 1
    return total


def compute_metrics(dna: str, obj: OptimizationObjectives, table: CodonTable) -> SequenceMetrics:
    lo, hi = local_gc_extremes(dna, obj.gc_window)
    return SequenceMetrics(
        length_nt=len(normalize_dna(dna)),
        gc_content=round(gc_content(dna), 4),
        cai=round(table.cai(dna), 4),
        max_homopolymer=max_homopolymer_run(dna),
        min_local_gc=round(lo, 4),
        max_local_gc=round(hi, 4),
        forbidden_motif_hits=_count_motifs(dna, obj.forbidden_motifs),
        restriction_site_hits=_count_motifs(dna, obj.avoid_restriction_sites),
    )


def penalty(metrics: SequenceMetrics, obj: OptimizationObjectives) -> float:
    """Lower is better."""
    codon_pen = obj.weight_codon * (1.0 - metrics.cai)

    gc_pen = 0.0
    g = metrics.gc_content
    if g < obj.target_gc_min:
        gc_pen += (obj.target_gc_min - g)
    elif g > obj.target_gc_max:
        gc_pen += (g - obj.target_gc_max)
    # local windows out of band
    if metrics.min_local_gc < obj.target_gc_min:
        gc_pen += (obj.target_gc_min - metrics.min_local_gc)
    if metrics.max_local_gc > obj.target_gc_max:
        gc_pen += (metrics.max_local_gc - obj.target_gc_max)
    gc_pen *= obj.weight_gc

    hp_pen = obj.weight_homopolymer * max(0, metrics.max_homopolymer - obj.max_homopolymer)
    motif_pen = obj.weight_motif * (metrics.forbidden_motif_hits + metrics.restriction_site_hits)
    return round(codon_pen + gc_pen + hp_pen + motif_pen, 6)


def _codons(protein: str) -> list[str]:
    return list(protein.upper())


def _baseline(protein: str, table: CodonTable) -> list[str]:
    return [table.best_codon(aa) for aa in _codons(protein)]


def _sampled(protein: str, table: CodonTable, rng: random.Random) -> list[str]:
    out: list[str] = []
    for aa in _codons(protein):
        codons = SYNONYMS[aa]
        weights = [max(table.fraction.get(c, 0.0), 1e-6) for c in codons]
        out.append(rng.choices(codons, weights=weights, k=1)[0])
    return out


def _stop_codon(table: CodonTable) -> str:
    return table.best_codon("*")


def _repair(
    codons: list[str],
    protein: str,
    obj: OptimizationObjectives,
    table: CodonTable,
    rng: random.Random,
) -> list[str]:
    """Iteratively swap synonymous codons to reduce penalty. Deterministic
    given rng. Protein identity is preserved because only synonymous codons
    are ever substituted."""
    aas = _codons(protein)

    def seq_of(cs: list[str]) -> str:
        return "".join(cs)

    for _ in range(6):  # bounded passes
        improved = False
        # evaluate current penalty
        cur_seq = seq_of(codons) + _stop_codon(table)
        cur_pen = penalty(compute_metrics(cur_seq, obj, table), obj)
        # try a synonymous swap at each position, accept if it lowers penalty
        order = list(range(len(codons)))
        rng.shuffle(order)
        for pos in order:
            aa = aas[pos]
            alts = [c for c in table.ranked_codons(aa) if c != codons[pos]]
            for alt in alts:
                trial = codons.copy()
                trial[pos] = alt
                trial_seq = seq_of(trial) + _stop_codon(table)
                trial_pen = penalty(compute_metrics(trial_seq, obj, table), obj)
                if trial_pen < cur_pen - 1e-9:
                    codons = trial
                    cur_pen = trial_pen
                    improved = True
                    break
        if not improved:
            break
    return codons


def _build_candidate(
    cid: str,
    codons: list[str],
    protein: str,
    obj: OptimizationObjectives,
    table: CodonTable,
    is_baseline: bool,
    note: str,
) -> OptimizationCandidate:
    cds_dna = "".join(codons) + _stop_codon(table)
    # safety: verify synonymous assembly really encodes the protein
    assert translate(cds_dna, to_stop=True) == protein.upper(), "optimizer broke protein identity"
    metrics = compute_metrics(cds_dna, obj, table)
    cds_rna = normalize_dna(cds_dna).replace("T", "U")
    return OptimizationCandidate(
        candidate_id=cid,
        cds_rna=cds_rna,
        metrics=metrics,
        objective_score=penalty(metrics, obj),
        is_baseline=is_baseline,
        note=note,
    )


def optimize(
    protein: str,
    objectives: OptimizationObjectives | None = None,
    table: CodonTable | None = None,
) -> OptimizationResult:
    """Optimize a protein (amino acids, NO stop) into a CDS including a stop."""
    obj = objectives or OptimizationObjectives()
    tbl = table or load_table(obj.host)
    protein = protein.upper()

    candidates: list[OptimizationCandidate] = []

    # baseline (max codon usage), then repaired baseline
    base = _baseline(protein, tbl)
    candidates.append(
        _build_candidate(
            "baseline_max_usage", base, protein, obj, tbl, True,
            "Highest-usage codon per residue (reference point).",
        )
    )
    base_rep = _repair(base.copy(), protein, obj, tbl, random.Random(obj.seed))
    candidates.append(
        _build_candidate(
            "baseline_repaired", base_rep, protein, obj, tbl, False,
            "Max-usage baseline after constraint-repair pass.",
        )
    )

    # seeded sampled + repaired alternatives
    n_alt = max(0, obj.n_candidates - 1)
    for i in range(n_alt):
        rng = random.Random(obj.seed + 1 + i)
        samp = _sampled(protein, tbl, rng)
        samp = _repair(samp, protein, obj, tbl, rng)
        candidates.append(
            _build_candidate(
                f"sampled_{i+1}", samp, protein, obj, tbl, False,
                f"Seeded weighted-sampling candidate #{i+1} (seed {obj.seed + 1 + i}), repaired.",
            )
        )

    # selection: minimum penalty; tie-break by higher CAI then candidate_id
    best = min(
        candidates,
        key=lambda c: (c.objective_score, -c.metrics.cai, c.candidate_id),
    )
    baseline = next(c for c in candidates if c.is_baseline)
    rationale = (
        f"Selected {best.candidate_id} with penalty {best.objective_score} "
        f"(CAI {best.metrics.cai}, GC {best.metrics.gc_content}, "
        f"max homopolymer {best.metrics.max_homopolymer}). "
        f"Baseline {baseline.candidate_id} penalty {baseline.objective_score} "
        f"(CAI {baseline.metrics.cai}, GC {baseline.metrics.gc_content}). "
        "Lower penalty is better; see candidate comparison for measured tradeoffs."
    )

    return OptimizationResult(
        objectives=obj,
        candidates=candidates,
        selected_id=best.candidate_id,
        selection_rationale=rationale,
        optimizer_version=OPTIMIZER_VERSION,
        codon_table_version=f"{tbl.key}:{tbl.version}",
        seed=obj.seed,
    )
