"""Load bundled codon-usage snapshots and compute CAI weights."""
from __future__ import annotations

import json
import pathlib
from functools import lru_cache

from .genetic_code import CODON_TABLE, SYNONYMS

_DATA = pathlib.Path(__file__).parent / "data" / "codon_tables"


class CodonTable:
    def __init__(self, raw: dict):
        self.key: str = raw["key"]
        self.name: str = raw["name"]
        self.source: str = raw["source"]
        self.version: str = raw["version"]
        self.note: str = raw["note"]
        # fraction per codon within its amino-acid family
        self.fraction: dict[str, float] = {}
        for aa, fam in raw["fractions_by_aa"].items():
            for codon, frac in fam.items():
                self.fraction[codon] = float(frac)
        # relative adaptiveness w = fraction / max(fraction in family); used for CAI
        self.w: dict[str, float] = {}
        for aa, codons in SYNONYMS.items():
            fr = {c: self.fraction.get(c, 0.0) for c in codons}
            mx = max(fr.values()) if fr else 0.0
            for c in codons:
                self.w[c] = (fr[c] / mx) if mx > 0 else 0.0

    def best_codon(self, aa: str) -> str:
        return max(SYNONYMS[aa], key=lambda c: self.fraction.get(c, 0.0))

    def ranked_codons(self, aa: str) -> list[str]:
        return sorted(SYNONYMS[aa], key=lambda c: self.fraction.get(c, 0.0), reverse=True)

    def cai(self, dna: str) -> float:
        """Codon Adaptation Index (geometric mean of relative adaptiveness).
        Stop codons and single-codon families (M, W) are excluded, per the
        standard CAI definition."""
        from .genetic_code import normalize_dna

        s = normalize_dna(dna)
        import math

        logs: list[float] = []
        for i in range(0, len(s) - len(s) % 3, 3):
            codon = s[i : i + 3]
            aa = CODON_TABLE.get(codon)
            if aa in (None, "*"):
                continue
            if len(SYNONYMS[aa]) <= 1:  # M, W carry no information
                continue
            wv = self.w.get(codon, 0.0)
            if wv <= 0:
                wv = 1e-6
            logs.append(math.log(wv))
        if not logs:
            return 0.0
        return math.exp(sum(logs) / len(logs))


@lru_cache(maxsize=None)
def load_table(key: str) -> CodonTable:
    path = _DATA / f"{key}.json"
    if not path.exists():
        raise KeyError(f"unknown codon table {key!r}; available: {available()}")
    return CodonTable(json.loads(path.read_text()))


def available() -> list[str]:
    return sorted(p.stem for p in _DATA.glob("*.json"))


def load_table_from_dict(raw: dict) -> CodonTable:
    """For user-supplied tables (same JSON schema as bundled files)."""
    return CodonTable(raw)
