"""Typed schemas shared across module boundaries.

These Pydantic models are the *contract* between the pipeline stages
(interpret -> references -> literature -> assemble -> transcript ->
optimize -> validate -> export). Every stage consumes and returns
these objects so that provenance and status travel with the data.

Design rule enforced structurally here:
  * The *intended protein sequence*, the *optimized CDS*, and the
    *full mature RNA* are distinct fields. They are never the same
    string reused, so one can never be silently mistaken for another.
"""
from __future__ import annotations

from datetime import date
from enum import Enum
from typing import Literal, Optional

from pydantic import BaseModel, Field


# --------------------------------------------------------------------------
# Provenance primitives
# --------------------------------------------------------------------------
class Provenance(str, Enum):
    """Where a piece of information came from. Central to traceability."""

    USER = "user"  # explicit user choice / supplied value
    REFERENCE = "reference"  # fact extracted from an authoritative sequence DB
    LITERATURE = "literature"  # fact extracted from a primary source
    ALGORITHM = "algorithm"  # output of a documented deterministic method
    PREDICTION = "prediction"  # computational prediction (not measured)
    ASSUMPTION = "assumption"  # a clearly-labeled design assumption
    FIXTURE = "fixture"  # bundled demo data, NOT live research


class DesignStatus(str, Enum):
    PROVISIONAL = "provisional"  # has unresolved dependencies / assumptions
    VALIDATED = "computationally_validated"  # passed all hard checks
    BLOCKED = "blocked"  # cannot proceed: essential input unresolved


class Assumption(BaseModel):
    """An explicitly recorded design assumption."""

    field: str
    value: str
    rationale: str
    provenance: Provenance = Provenance.ASSUMPTION


# --------------------------------------------------------------------------
# Request
# --------------------------------------------------------------------------
class ConstructContext(BaseModel):
    """Optional structured context. Encoded-protein species is kept
    separate from the expression host on purpose."""

    protein_species: Optional[str] = None  # species of the ENCODED protein
    expression_host: Optional[str] = None  # organism whose codon usage to use
    cell_type: Optional[str] = None
    in_vivo: Optional[bool] = None  # True=in vivo, False=in vitro, None=unspecified
    delivery: Optional[str] = None  # e.g. "LNP", "electroporation"
    localization: Optional[str] = None  # e.g. "secreted", "cytoplasmic"
    max_total_nt: Optional[int] = None  # user size constraint (biological)


class ConstructRequest(BaseModel):
    text: str = ""  # natural-language description
    protein_identifier: Optional[str] = None  # name or accession
    pasted_protein_sequence: Optional[str] = None  # amino acids
    user_nucleotide_sequence: Optional[str] = None  # optional user CDS
    preferred_elements: dict[str, str] = Field(default_factory=dict)
    context: ConstructContext = Field(default_factory=ConstructContext)
    # If the user edited the interpreted architecture directly, it rides here:
    architecture_override: Optional["Architecture"] = None


# --------------------------------------------------------------------------
# Clarification
# --------------------------------------------------------------------------
class ClarificationOption(BaseModel):
    label: str
    value: str


class ClarificationQuestion(BaseModel):
    key: str
    question: str
    options: list[ClarificationOption] = Field(default_factory=list)
    impact: str  # why this question matters to the design


class ClarificationRound(BaseModel):
    """At most two questions, exactly one round."""

    questions: list[ClarificationQuestion] = Field(default_factory=list, max_length=2)


# --------------------------------------------------------------------------
# Architecture (structured construct spec, BEFORE nucleotides)
# --------------------------------------------------------------------------
ComponentKind = Literal[
    "protein", "domain", "signal_peptide", "tag", "linker", "localization"
]


class Component(BaseModel):
    kind: ComponentKind
    name: str
    sequence: str = ""  # amino acids for this component (may be "" if unresolved)
    species: Optional[str] = None
    accession: Optional[str] = None
    accession_version: Optional[str] = None
    residue_start: Optional[int] = None  # 1-based inclusive, within source protein
    residue_end: Optional[int] = None
    placement: Optional[str] = None  # e.g. "N", "C", "internal"
    provenance: Provenance = Provenance.ASSUMPTION
    source_id: Optional[str] = None  # link into evidence records
    resolved: bool = False  # False => sequence unknown / provisional
    note: str = ""


class Substitution(BaseModel):
    """User-requested point substitution or deletion, 1-based."""

    position: int
    from_aa: str
    to_aa: str  # "" denotes deletion


class Architecture(BaseModel):
    title: str = "Untitled construct"
    components: list[Component] = Field(default_factory=list)
    fusion_order: list[int] = Field(default_factory=list)  # indices into components
    substitutions: list[Substitution] = Field(default_factory=list)
    is_precursor: bool = False  # contains signal peptide / propeptide
    mature_note: str = ""  # how mature chain differs from precursor
    unresolved: list[str] = Field(default_factory=list)  # human-readable blockers
    assumptions: list[Assumption] = Field(default_factory=list)
    supported: bool = True  # False => unsupported architecture, explained in notes

    def assembled_protein(self) -> str:
        """Concatenate component sequences in fusion order (no substitutions)."""
        order = self.fusion_order or list(range(len(self.components)))
        return "".join(self.components[i].sequence for i in order)


# --------------------------------------------------------------------------
# Evidence
# --------------------------------------------------------------------------
class EvidenceRecord(BaseModel):
    source_id: str
    title: str
    url: str = ""
    doi: Optional[str] = None
    pmid: Optional[str] = None
    accession: Optional[str] = None
    retrieved: Optional[date] = None
    finding: str = ""  # relevant finding (a fact from the source)
    experimental_context: str = ""  # species/cell the source actually used
    exact_sequence_available: bool = False
    supports_decision: str = ""  # which design decision this backs
    match_quality: Literal["exact", "analogous", "background", "none"] = "background"
    provenance: Provenance = Provenance.LITERATURE
    conflicts_with: list[str] = Field(default_factory=list)  # other source_ids


class RetrievalState(BaseModel):
    """Preserves the success/failure of each retrieval integration."""

    integration: str
    attempted: bool = False
    live: bool = False  # True if a network call actually happened
    ok: bool = False
    detail: str = ""


# --------------------------------------------------------------------------
# Transcript elements + transcript
# --------------------------------------------------------------------------
class CapSpec(BaseModel):
    name: str = "Cap1"  # metadata only; NOT inserted into sequence
    chemistry: str = ""
    provenance: Provenance = Provenance.USER
    source_id: Optional[str] = None
    note: str = ""


class PolyASpec(BaseModel):
    mode: Literal["encoded", "post_transcriptional"] = "encoded"
    length: int = 120
    provenance: Provenance = Provenance.ASSUMPTION
    source_id: Optional[str] = None
    note: str = ""


class TranscriptElement(BaseModel):
    """A versioned, source-backed UTR (or similar) from the library."""

    element_id: str
    role: Literal["5utr", "3utr", "kozak"]
    name: str
    sequence: str  # RNA or DNA letters (A/C/G/T or U)
    version: str
    source_id: str
    relevance: str
    uncertainty: str  # explicit statement of how well-supported this is
    provenance: Provenance = Provenance.REFERENCE


class ModNucleoside(BaseModel):
    name: str  # e.g. "N1-methylpseudouridine"
    replaces: str  # e.g. "U"
    provenance: Provenance = Provenance.USER
    note: str = ""


class Transcript(BaseModel):
    """Distinct sequence objects, never conflated."""

    intended_protein: str  # amino acids the construct must encode
    optimized_cds: str  # RNA CDS INCLUDING stop codon (letters A/C/G/U)
    five_utr: TranscriptElement | None = None
    kozak: TranscriptElement | None = None
    three_utr: TranscriptElement | None = None
    polyA: PolyASpec = Field(default_factory=PolyASpec)
    cap: CapSpec = Field(default_factory=CapSpec)
    mod_nucleosides: list[ModNucleoside] = Field(default_factory=list)

    def mature_rna(self) -> str:
        """Full exported mature RNA (no promoter, no cap letters).
        poly(A) letters are included ONLY when mode == 'encoded'."""
        parts: list[str] = []
        if self.five_utr:
            parts.append(_to_rna(self.five_utr.sequence))
        if self.kozak:
            parts.append(_to_rna(self.kozak.sequence))
        parts.append(_to_rna(self.optimized_cds))
        if self.three_utr:
            parts.append(_to_rna(self.three_utr.sequence))
        if self.polyA.mode == "encoded":
            parts.append("A" * self.polyA.length)
        return "".join(parts)


def _to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


# --------------------------------------------------------------------------
# Optimization
# --------------------------------------------------------------------------
class OptimizationObjectives(BaseModel):
    host: str = "h_sapiens"  # codon table key
    target_gc_min: float = 0.40
    target_gc_max: float = 0.70
    gc_window: int = 40  # local GC window (nt)
    max_homopolymer: int = 9  # longest allowed run of identical nt
    forbidden_motifs: list[str] = Field(default_factory=list)
    avoid_restriction_sites: list[str] = Field(default_factory=list)
    weight_codon: float = 1.0
    weight_gc: float = 1.0
    weight_homopolymer: float = 1.0
    weight_motif: float = 1.0
    n_candidates: int = 3
    seed: int = 12345


class SequenceMetrics(BaseModel):
    length_nt: int
    gc_content: float
    cai: float  # Codon Adaptation Index vs host table
    max_homopolymer: int
    min_local_gc: float
    max_local_gc: float
    forbidden_motif_hits: int
    restriction_site_hits: int


class OptimizationCandidate(BaseModel):
    candidate_id: str
    cds_rna: str  # includes stop codon
    metrics: SequenceMetrics
    objective_score: float
    is_baseline: bool = False
    note: str = ""


class OptimizationResult(BaseModel):
    objectives: OptimizationObjectives
    candidates: list[OptimizationCandidate]
    selected_id: str
    selection_rationale: str
    optimizer_version: str
    codon_table_version: str
    seed: int


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
class ValidationItem(BaseModel):
    check: str
    severity: Literal["error", "warning", "unevaluated", "pass"]
    message: str


class SizeReport(BaseModel):
    protein_len_aa: int
    cds_len_nt: int
    cds_includes_stop: bool
    five_utr_len_nt: int
    three_utr_len_nt: int
    polyA_encoded_len_nt: int  # 0 if post-transcriptional
    polyA_mode: str
    total_rna_len_nt: int
    software_input_limit_nt: int
    biological_guidance: str
    guidance_source_id: Optional[str] = None
    guidance_is_uncertain: bool = True
    component_contributions: dict[str, int] = Field(default_factory=dict)


class ValidationReport(BaseModel):
    items: list[ValidationItem] = Field(default_factory=list)
    translated_protein: str = ""
    size: SizeReport | None = None

    @property
    def errors(self) -> list[ValidationItem]:
        return [i for i in self.items if i.severity == "error"]

    @property
    def warnings(self) -> list[ValidationItem]:
        return [i for i in self.items if i.severity == "warning"]

    @property
    def unevaluated(self) -> list[ValidationItem]:
        return [i for i in self.items if i.severity == "unevaluated"]

    @property
    def passed(self) -> bool:
        return len(self.errors) == 0


# --------------------------------------------------------------------------
# Top-level design object + manifest
# --------------------------------------------------------------------------
class Design(BaseModel):
    design_id: str
    request: ConstructRequest
    architecture: Architecture
    evidence: list[EvidenceRecord] = Field(default_factory=list)
    retrieval_states: list[RetrievalState] = Field(default_factory=list)
    transcript: Transcript | None = None
    optimization: OptimizationResult | None = None
    validation: ValidationReport | None = None
    status: DesignStatus = DesignStatus.PROVISIONAL
    created: date = Field(default_factory=date.today)
    is_demo_fixture: bool = False  # True => uses bundled data, not live research


# Resolve forward reference ConstructRequest.architecture_override
ConstructRequest.model_rebuild()
