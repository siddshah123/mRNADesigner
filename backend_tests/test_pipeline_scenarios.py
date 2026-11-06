"""Scenario tests covering the required delivery matrix.

Single reference protein · tagged protein · two-domain fusion w/ linker ·
ambiguous identity · unresolved isoform · length-limit violation ·
retrieval failure · translation mismatch · export + re-import.
"""
import io

import pytest
from Bio import SeqIO

from mrna_designer import export, literature
from mrna_designer import references as R
from mrna_designer.genetic_code import translate
from mrna_designer.interpret import interpret
from mrna_designer.pipeline import run_pipeline
from mrna_designer.schemas import (
    ConstructContext,
    ConstructRequest,
    DesignStatus,
    OptimizationObjectives,
)


def _req(**kw):
    return ConstructRequest(**kw)


# ---- 1. single reference protein -----------------------------------------
def test_single_reference_protein():
    d, clar = run_pipeline(_req(protein_identifier="demo_reporter"))
    assert d.transcript is not None
    assert d.validation.passed  # no hard errors
    assert translate(d.transcript.optimized_cds, to_stop=True) == d.transcript.intended_protein
    assert d.is_demo_fixture is True
    assert d.status == DesignStatus.PROVISIONAL  # demo => provisional, not "validated"


# ---- 2. tagged protein ----------------------------------------------------
def test_tagged_protein_appends_tag():
    d, _ = run_pipeline(_req(protein_identifier="demo_reporter",
                             text="demo_reporter with a C-terminal His tag"))
    kinds = [c.kind for c in d.architecture.components]
    assert "tag" in kinds
    tag = next(c for c in d.architecture.components if c.kind == "tag")
    assert tag.sequence == "HHHHHH"
    # tag residues must appear at the end of the intended protein (C-terminal)
    assert d.transcript.intended_protein.endswith("HHHHHH")


# ---- 3. two-domain fusion with a linker -----------------------------------
def test_two_domain_fusion_with_linker():
    d, _ = run_pipeline(_req(
        protein_identifier="demo_reporter fused to demo_secreted_precursor",
        text="fusion with a linker"))
    kinds = [c.kind for c in d.architecture.components]
    assert "linker" in kinds
    assert kinds.count("protein") >= 2
    # junction validity is checked
    assert any(i.check == "fusion_junctions" and i.severity == "pass"
               for i in d.validation.items)
    # linker sequence is present inside the protein
    linker = next(c for c in d.architecture.components if c.kind == "linker")
    assert linker.sequence in d.transcript.intended_protein


# ---- 4. ambiguous protein identity ---------------------------------------
def test_ambiguous_identity_triggers_clarification():
    d, clar = run_pipeline(_req(protein_identifier="interleukin-2"))
    assert d.status == DesignStatus.BLOCKED
    assert d.transcript is None  # no CDS invented
    assert any(q.key == "protein_identity" for q in clar.questions)
    assert d.architecture.unresolved


# ---- 5. unresolved isoform / accession not retrievable offline -----------
def test_unresolved_accession_stays_provisional_no_invention():
    # a valid-looking UniProt accession, but network off => cannot resolve
    d, clar = run_pipeline(_req(protein_identifier="P60568"), allow_network=False)
    assert d.status == DesignStatus.BLOCKED
    assert d.transcript is None
    comp = d.architecture.components[0]
    assert comp.resolved is False
    assert comp.sequence == ""  # nothing fabricated


# ---- 6. length-limit violation -------------------------------------------
def test_biological_size_warning_is_soft_and_uncertain():
    # a big user 5'UTR pushes total over the LNP soft ceiling (~5000 nt)
    d, _ = run_pipeline(
        _req(protein_identifier="demo_reporter",
             preferred_elements={"five_utr": "A" * 6000}),
        size_profile_id="lnp_ivt_mrna")
    warn = [i for i in d.validation.items if i.check == "biological_size_guidance"]
    assert warn and warn[0].severity == "warning"
    assert "uncertain" in warn[0].message.lower()
    # it is NOT a hard error
    assert all(i.check != "biological_size_guidance" or i.severity != "error"
               for i in d.validation.items)


def test_software_size_limit_is_hard_error():
    d, _ = run_pipeline(
        _req(protein_identifier="demo_reporter",
             preferred_elements={"five_utr": "A" * 60000}))
    hard = [i for i in d.validation.items if i.check == "software_size_limit"]
    assert hard and hard[0].severity == "error"
    assert "contribution" in hard[0].message.lower() or "5utr=" in hard[0].message


# ---- 7. retrieval failure -------------------------------------------------
def test_uniprot_retrieval_failure_is_recorded_not_faked():
    ref = R.fetch_uniprot("P60568", allow_network=True)
    assert ref.retrieval.attempted
    # in this environment the proxy blocks it; either way, no sequence is faked on failure
    if not ref.retrieval.ok:
        assert ref.sequence == ""
        assert ref.retrieval.detail


def test_uniprot_name_search_resolves_when_network_mocked(monkeypatch):
    """With network ON and a reachable UniProt, a plain NAME resolves to a real
    sequence (deterministic research, no LLM, no fabrication)."""
    import httpx

    fake_fasta = (">sp|P60568|IL2_HUMAN Interleukin-2 OS=Homo sapiens OX=9606 "
                  "GN=IL2 PE=1 SV=1\nMAPTSSSTKK\nTQLQLEHLLL\n")

    class _Resp:
        status_code = 200
        text = fake_fasta

    def fake_get(url, **kw):
        assert "rest.uniprot.org" in url
        return _Resp()

    monkeypatch.setattr(httpx, "get", fake_get)

    ref = R.search_uniprot("interleukin-2", species="Homo sapiens", allow_network=True)
    assert ref.retrieval.ok and ref.accession == "P60568"
    assert ref.sequence == "MAPTSSSTKKTQLQLEHLLL"
    assert ref.provenance.value == "reference"

    # and the interpreter uses it to resolve a named protein end-to-end
    arch, _ = interpret(_req(protein_identifier="interleukin-2"), allow_network=True)
    comp = arch.components[0]
    assert comp.resolved and comp.sequence == "MAPTSSSTKKTQLQLEHLLL"
    assert comp.accession == "P60568"


def test_bare_name_offline_stays_unresolved_not_invented():
    arch, _ = interpret(_req(text="IL4 and lumazine synthase"), allow_network=False)
    prot = [c for c in arch.components if c.kind == "protein"]
    assert len(prot) == 2
    assert all(not c.resolved and c.sequence == "" for c in prot)  # nothing invented


def test_live_literature_unavailable_state_preserved():
    ev, state = literature.retrieve(["utr"], allow_network=True)
    assert state.attempted and state.ok is False and state.live is False
    assert "not enabled" in state.detail.lower()


# ---- 8. translation mismatch (hard failure) ------------------------------
def test_translation_mismatch_is_hard_error():
    from mrna_designer.assemble import assemble
    from mrna_designer.transcript import build_transcript
    from mrna_designer.validate import validate_design

    arch, _ = interpret(_req(protein_identifier="demo_reporter"))
    asm = assemble(arch)
    # deliberately corrupt the CDS so translation != intended protein
    bad_cds = "AUG" + "AAA" * 5 + "UAA"  # encodes MKKKKK*, not the reporter
    tr = build_transcript(asm.precursor_protein, bad_cds)
    rep = validate_design(intended_protein=asm.precursor_protein, assembly=asm,
                          transcript=tr, objectives=OptimizationObjectives())
    assert any(i.check == "translation_identity" and i.severity == "error"
               for i in rep.items)
    assert not rep.passed


def test_internal_stop_is_hard_error():
    from mrna_designer.transcript import build_transcript
    from mrna_designer.validate import validate_design

    prot = "MKKK"
    bad = "AUGUAAAAAAAAUAA"  # has an internal stop
    tr = build_transcript(prot, bad)
    rep = validate_design(intended_protein=prot, assembly=None, transcript=tr,
                          objectives=OptimizationObjectives())
    assert any(i.check == "internal_stop" and i.severity == "error" for i in rep.items)


# ---- 9. export + re-import round trip -------------------------------------
def test_genbank_export_roundtrips_with_features():
    d, _ = run_pipeline(_req(protein_identifier="demo_secreted_precursor",
                             context=ConstructContext(expression_host="h_sapiens")))
    rt = export.roundtrip_check(d)
    assert rt["genbank_roundtrip"] == "pass"
    assert rt["fasta_roundtrip"] == "pass"
    # re-parse and verify a CDS feature with correct translation exists
    gb = export.genbank(d)
    rec = next(SeqIO.parse(io.StringIO(gb), "genbank"))
    cds = [f for f in rec.features if f.type == "CDS"]
    assert cds
    assert cds[0].qualifiers["translation"][0] == d.transcript.intended_protein
    # sig_peptide feature present (precursor)
    assert any(f.type == "sig_peptide" for f in rec.features)


def test_no_promoter_and_cap_is_metadata_only():
    d, _ = run_pipeline(_req(protein_identifier="demo_reporter"))
    rna = d.transcript.mature_rna()
    # cap name must not appear as letters; sequence is pure A/C/G/U
    assert set(rna) <= set("ACGU")
    gb = export.genbank(d)
    rec = next(SeqIO.parse(io.StringIO(gb), "genbank"))
    # no promoter FEATURE in the mature-mRNA export (the comment may mention it)
    assert not any("promoter" in f.type.lower() for f in rec.features)
    assert "METADATA ONLY" in gb  # cap disclosed as metadata in the comment


def test_polyA_encoded_vs_length_consistency():
    d, _ = run_pipeline(_req(protein_identifier="demo_reporter"))
    s = d.validation.size
    rna = d.transcript.mature_rna()
    if d.transcript.polyA.mode == "encoded":
        assert rna.endswith("A" * s.polyA_encoded_len_nt)
        assert s.polyA_encoded_len_nt == d.transcript.polyA.length
