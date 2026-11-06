"""The standard genetic code (NCBI transl_table=1) and basic sequence ops.

Deterministic, dependency-free. Works in DNA letters internally (A/C/G/T);
RNA is produced only at export by T->U substitution.
"""
from __future__ import annotations

# codon (DNA) -> amino acid (1-letter), '*' = stop
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# amino acid -> list of synonymous codons (DNA)
SYNONYMS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    SYNONYMS.setdefault(_aa, []).append(_codon)

AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")
STOP = "*"


def normalize_dna(seq: str) -> str:
    return seq.strip().upper().replace("U", "T").replace(" ", "").replace("\n", "")


def translate(dna: str, to_stop: bool = False) -> str:
    """Translate DNA (or RNA) in frame 0. Returns amino acids with '*' for stop.
    Raises ValueError if length is not a multiple of 3 or an invalid codon."""
    s = normalize_dna(dna)
    if len(s) % 3 != 0:
        raise ValueError(f"sequence length {len(s)} is not a multiple of 3")
    aas: list[str] = []
    for i in range(0, len(s), 3):
        codon = s[i : i + 3]
        if codon not in CODON_TABLE:
            raise ValueError(f"invalid codon {codon!r} at nt {i}")
        aa = CODON_TABLE[codon]
        if aa == STOP and to_stop:
            break
        aas.append(aa)
    return "".join(aas)


def is_valid_protein(seq: str) -> bool:
    return len(seq) > 0 and all(c in AA_ALPHABET for c in seq.upper())


def gc_content(seq: str) -> float:
    s = normalize_dna(seq)
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s)


def max_homopolymer_run(seq: str) -> int:
    s = normalize_dna(seq)
    best = run = 0
    prev = ""
    for c in s:
        run = run + 1 if c == prev else 1
        prev = c
        best = max(best, run)
    return best


def local_gc_extremes(seq: str, window: int) -> tuple[float, float]:
    """Min and max GC fraction over sliding windows. Returns (min, max)."""
    s = normalize_dna(seq)
    if len(s) < window or window <= 0:
        g = gc_content(s)
        return (g, g)
    gc = [1 if c in "GC" else 0 for c in s]
    # prefix sums
    pref = [0]
    for v in gc:
        pref.append(pref[-1] + v)
    lo, hi = 1.0, 0.0
    for i in range(0, len(s) - window + 1):
        frac = (pref[i + window] - pref[i]) / window
        lo = min(lo, frac)
        hi = max(hi, frac)
    return (lo, hi)
