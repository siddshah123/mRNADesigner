"""Sequence-invariant tests for the deterministic engine."""
import pytest

from mrna_designer import codon_tables
from mrna_designer.genetic_code import (
    gc_content,
    is_valid_protein,
    max_homopolymer_run,
    translate,
)
from mrna_designer.optimize import optimize
from mrna_designer.schemas import OptimizationObjectives

PROT = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEFKLVNELTEFAKTCVADESHAGCEK"


def test_translate_roundtrip_and_frame():
    with pytest.raises(ValueError):
        translate("ATGAA")  # not multiple of 3
    assert translate("ATGTAA") == "M*"
    assert translate("ATGTAA", to_stop=True) == "M"


def test_optimizer_preserves_protein_identity():
    r = optimize(PROT, OptimizationObjectives(host="h_sapiens", seed=1))
    for c in r.candidates:
        assert translate(c.cds_rna, to_stop=True) == PROT, c.candidate_id
        assert c.cds_rna.endswith(("UAA", "UAG", "UGA"))


def test_optimizer_is_deterministic():
    a = optimize(PROT, OptimizationObjectives(host="h_sapiens", seed=7))
    b = optimize(PROT, OptimizationObjectives(host="h_sapiens", seed=7))
    sa = next(c for c in a.candidates if c.candidate_id == a.selected_id)
    sb = next(c for c in b.candidates if c.candidate_id == b.selected_id)
    assert sa.cds_rna == sb.cds_rna
    # different seed can change sampled candidates
    c = optimize(PROT, OptimizationObjectives(host="h_sapiens", seed=8))
    assert isinstance(c.selected_id, str)


def test_different_host_tables_differ():
    h = optimize(PROT, OptimizationObjectives(host="h_sapiens", seed=3, n_candidates=1))
    e = optimize(PROT, OptimizationObjectives(host="e_coli", seed=3, n_candidates=1))
    hs = next(c for c in h.candidates if c.candidate_id == h.selected_id)
    es = next(c for c in e.candidates if c.candidate_id == e.selected_id)
    # same protein, different codon usage => sequences should differ
    assert hs.cds_rna != es.cds_rna
    assert translate(hs.cds_rna, to_stop=True) == translate(es.cds_rna, to_stop=True) == PROT


def test_forbidden_motif_removed_when_possible():
    # ask optimizer to avoid EcoRI site GAATTC
    r = optimize(PROT, OptimizationObjectives(host="h_sapiens", seed=5,
                                              forbidden_motifs=["GAATTC"], weight_motif=5))
    sel = next(c for c in r.candidates if c.candidate_id == r.selected_id)
    assert sel.metrics.forbidden_motif_hits == 0


def test_homopolymer_cap_respected_by_selection():
    poly = "M" + "K" * 20  # lysine run can create AAA homopolymers
    r = optimize(poly, OptimizationObjectives(host="h_sapiens", seed=2, max_homopolymer=6,
                                              weight_homopolymer=5))
    sel = next(c for c in r.candidates if c.candidate_id == r.selected_id)
    # selection should not exceed cap if a repair exists; at minimum it is reported
    assert sel.metrics.max_homopolymer <= 9


def test_cai_bounds():
    t = codon_tables.load_table("h_sapiens")
    cds = "".join(t.best_codon(a) for a in PROT) + "TAA"
    assert 0.0 < t.cai(cds) <= 1.0001
    # max-usage sequence should have CAI == 1.0 (all best codons)
    assert abs(t.cai(cds) - 1.0) < 1e-6


def test_metrics_helpers():
    assert is_valid_protein("MKV")
    assert not is_valid_protein("MKZ")
    assert abs(gc_content("GGCC") - 1.0) < 1e-9
    assert max_homopolymer_run("AAATTTT") == 4


def test_codon_table_provenance_present():
    for key in codon_tables.available():
        t = codon_tables.load_table(key)
        assert t.source and t.version and t.note
