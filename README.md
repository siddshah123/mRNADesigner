# mRNADesigner

Convert a natural-language protein-construct request into a **traceable,
computationally-validated mRNA design**. Resolve the protein architecture,
attach supporting evidence, assemble the amino-acid sequence, perform
deterministic synonymous codon optimization, add transcript elements, validate
sequence invariants, and export annotated sequences with an evidence +
validation report.

> **Central principle — traceability.** Every amino acid, transcript element,
> optimization decision, and validation claim is attributable to a reference, an
> explicit user choice, a documented algorithm, or a clearly-labeled assumption.
> The app never invents a protein sequence, never fabricates a citation, and
> never presents a computational check as experimental proof.

This is the smallest scientifically useful version, built so the architecture
expands cleanly.

---

## Quick start

```bash
pip install -r requirements.txt
# (one-time) regenerate bundled codon tables — already committed, optional:
python -m mrna_designer._build_codon_tables

./run.sh                 # or: python -m uvicorn mrna_designer.api:app --reload
# open http://127.0.0.1:8000
```

Run the tests:

```bash
python -m pytest backend_tests/ -q
```

---

## Architecture

**Stack (local-first).** Python sequence-processing backend + a dependency-free
vanilla-JS browser UI, served by one FastAPI process. Pydantic v2 provides typed
schemas at every module boundary (`mrna_designer/schemas.py`). BioPython writes
GenBank and performs the export round-trip parse. **The language model is not in
the sequence path** — it is an optional hook for interpreting free text and
explaining evidence; all nucleotide generation is deterministic Python.

**Separated modules** (each independently testable):

| Stage | Module | Responsibility |
|-------|--------|----------------|
| Request interpretation | `interpret.py` | free text / id / pasted seq → structured `Architecture`; ≤2 clarification questions |
| Reference retrieval | `references.py` | bundled tags/linkers/demo proteins + **live UniProt fetch** with recorded failure state |
| Literature retrieval | `literature.py` | bounded topic-matched evidence; failure/unavailable state preserved |
| Construct assembly | `assemble.py` | components → precursor/mature protein, junction spans, substitutions |
| Transcript | `transcript.py` | versioned, source-backed UTR/Kozak/cap/poly(A) library + transcript build |
| Optimization | `optimize.py` | deterministic, seeded, multi-objective codon optimization |
| Validation | `validate.py` | hard failures + explicit size report (software limit vs biological guidance) |
| Export | `export.py` | RNA/CDS/protein FASTA, GenBank, JSON manifest, report + round-trip check |
| Orchestration | `pipeline.py` | wires the stages, sets final status |
| HTTP API + UI | `api.py`, `frontend/` | FastAPI routes + single-page interface |

Data flows as immutable typed objects that keep the **intended protein**, the
**optimized CDS**, and the **full mature RNA** as distinct fields so one is never
silently mistaken for another.

---

## Scientific assumptions (made explicit)

1. **Encoded-protein species ≠ expression host.** Two independent fields. A human
   protein is *not* assumed to be optimized for human cells; the host sets the
   codon table, the protein's species sets identity.
2. **Codon tables are bundled snapshots** (`data/codon_tables/*.json`): approximate
   relative synonymous codon usage consistent with the Kazusa database, each
   carrying `source`, `version`, and an explicit *"verify before production"*
   caveat. Upload your own table (same JSON schema) to override.
3. **Transcript elements are a small, versioned, source-backed library**
   (`data/references/transcript_elements.json`). The Kozak consensus is from
   Kozak (1987). The bundled UTRs are **synthetic placeholders**, labeled as such
   with an uncertainty statement — there is no "universally optimal" UTR or mRNA
   size. Supply your own UTRs for real work.
4. **poly(A): encoded vs post-transcriptional** are distinct. Encoded `A`s appear
   in the exported sequence; a post-transcriptional tail is recorded as intended
   length only, and the reported length matches the mode.
5. **Cap chemistry and modified nucleosides are metadata only** — never letters in
   FASTA, and disclosed as metadata in the GenBank comment.
6. **Retrieved content is data, never instructions.** Live UniProt/NCBI are
   contacted only when you opt in (`allow_network`); failures are recorded, not
   faked. Bundled demo proteins are **synthetic**, flagged everywhere, and are not
   real database entries.

### Status semantics

- **computationally_validated** — all hard checks pass, no unresolved inputs or
  assumptions, not demo data.
- **provisional** — assumptions recorded, demo data used, or non-fatal issues.
- **blocked** — an essential identity/sequence is unresolved; **no CDS is invented**.

---

## Deterministic codon optimization

Documented multi-objective strategy (`optimize.py`), not a black box:

1. **Candidates** — a *baseline* (highest-usage codon per residue, i.e. max codon
   adaptation — the reference point, *not* assumed best), plus seeded
   weighted-sampling alternatives (`Random(seed+i)`).
2. **Repair pass** on every candidate, only ever swapping a codon for a
   **synonymous** one (so the encoded protein is invariant by construction):
   breaks homopolymer runs, removes forbidden motifs / restriction sites, pulls
   local-GC windows into band.
3. **Scoring** — weighted penalty (lower = better) over codon usage `(1 − CAI)`,
   global + local GC deviation, homopolymer excess, motif/site hits.
4. **Selection** — minimum penalty; the rationale and the full candidate
   comparison (measured tradeoffs) are reported.

**Reproducible:** identical `(protein, objectives, seed)` reproduce byte-identical
output. Optimizer version, codon-table version, and seed are saved in the result.

Configurable objectives are visible in the UI: host codon table, global/local GC
band, homopolymer cap, forbidden motifs, candidate count, seed. Unsupported
analyses (e.g. RNA secondary structure — no folding tool is wired in this build)
are reported as **"unevaluated"**, never faked as a pass.

---

## Validation (hard failures)

Invalid alphabet · frameshift / incomplete triplets · internal stop codons ·
translation ≠ specified protein · incorrect fusion junctions · invalid feature
coordinates · explicit hard-constraint violations (homopolymer cap, forbidden
motifs) · software input-size limit.

The **size report** separates the *software input limit* (50 kb, a processing
ceiling) from *biological/delivery guidance* (a configurable, source-backed
profile, labeled **uncertain**, surfaced as a soft **warning** — never a hard
failure). When a construct exceeds a limit, the per-component nt contribution is
shown; nothing is silently truncated.

Every computational result carries the disclaimer that it is **not** evidence of
expression, secretion, activity, stability, or delivery.

---

## Exports

RNA FASTA · CDS FASTA · protein FASTA · annotated **GenBank** · JSON design
manifest · human-readable report. GenBank features annotate transcript elements
(5′UTR, Kozak, CDS, 3′UTR, encoded poly(A)) and protein components
(signal peptide, domains, linker, tag) with provenance.

- **No transcription promoter** in the mature-mRNA export. An optional IVT **DNA
  template** export is a *separately-labeled* format.
- RNA vs DNA representation is explicit (GenBank stores DNA letters with
  `molecule_type=mRNA` and a note).
- **MacVector:** GenBank is emitted as a standard `.gb` file, which MacVector
  imports. We do **not** rename another format with a MacVector extension; no
  documented native MacVector binary writer was used.
- Every export **round-trips** through BioPython (sequence identity + CDS feature
  coordinates); the result is shown in the UI and the `/api/design` response.

---

## Example workflow

Request: *"demo_reporter fused to demo_secreted_precursor with a C-terminal His
tag and a linker"*, expression host **h_sapiens**, forbidden motif `GAATTC`, LNP
size profile.

The pipeline splits the secreted precursor into a **signal peptide + mature
chain**, inserts a flexible `(GGGGS)×3` linker (cited: Chen, Zaro & Shen 2013),
appends a His6 tag, optimizes to a candidate with no `GAATTC` sites, assembles the
transcript (Kozak + synthetic UTRs + encoded poly(A)=120), validates (translation
identity ✓, junctions ✓), and exports a round-trip-verified GenBank.

A full rendered report is committed at [`docs/example_report.txt`](docs/example_report.txt).

---

## Integrations, optional dependencies, limitations

**Functioning now (offline):** request interpretation, architecture assembly
(incl. precursor/mature, fusions, tags, linkers, substitutions), deterministic
codon optimization, transcript assembly, full validation, all exports with
round-trip verification, bundled literature evidence, the complete browser UI.

**Only the description is required.** Every other field is optional; the more you
give, the more constrained the design. With less, the app researches what it can
and records explicit, labeled assumptions rather than guessing silently.

**Opt-in network (`allow_network=true`):** live **UniProt** retrieval — both an
accession fetch (`P60568`) and a **name search** (`IL4`, `lumazine synthase`),
which resolves a plain protein name to the top reviewed Swiss-Prot sequence
(deterministic lookup, **no LLM**). In a sandbox that blocks `rest.uniprot.org`
the failure is recorded honestly and the component is marked unresolved (design
stays provisional/blocked). No unpublished user sequence is sent anywhere unless
network retrieval is enabled; the response discloses every integration contacted
(`retrieval_states`).

**No LLM is connected.** The app requires no API key and no account access. The
language model is an optional interpretation hook only; all sequence resolution
is either a database lookup, a user-pasted sequence, or a labeled unresolved
dependency — never model-generated.

**Not wired in this build (shown as "unevaluated" / "unavailable", never faked):**
- Live literature search (NCBI E-utilities / Europe PMC) — the hook exists in
  `literature.py`; it currently returns bundled fixtures and records the
  unavailable state.
- RNA secondary-structure metrics (no folding engine bundled).
- Native MacVector binary export (GenBank import is the supported path).

**Known limitations:** bundled codon tables are approximate snapshots (verify
against Kazusa for production); bundled UTRs are synthetic placeholders; the demo
proteins are synthetic and benign; the free-text parser is a light rule-based
heuristic (an LLM interpreter can be slotted behind `interpret()` without
changing the deterministic sequence path); mid-component deletions re-derive
feature spans conservatively.

**Privacy:** API keys stay server-side (none are required offline). Fixtures are
visibly distinguished from live research throughout (status badge + report
banner + per-component provenance).

---

## Tests

`python -m pytest backend_tests/ -q` → **23 passed**. Coverage includes the
required matrix: single reference protein, tagged protein, two-domain fusion with
a linker, ambiguous identity, unresolved accession/isoform, length-limit
violation (soft biological warning + hard software error), retrieval failure,
translation mismatch, internal stop, and export + re-import round-trip — testing
sequence invariants and failure behavior, not just UI rendering.
