"""Builder that emits bundled codon-usage snapshots as JSON.

Run:  python -m mrna_designer._build_codon_tables

The numbers below are *approximate, widely-published* relative synonymous
codon usage fractions (fraction within each amino-acid family). They are
intended as a reproducible bundled default and carry an explicit caveat in
their `source`/`note` fields. For production work, replace with a freshly
downloaded table (Kazusa / HIVE-CUT) via the same JSON schema, or upload one
through the UI. Nothing in the app treats these as authoritative beyond the
stated provenance.
"""
import json
import pathlib

OUT = pathlib.Path(__file__).parent / "data" / "codon_tables"

# fraction within amino-acid family (values per family sum ~1.0)
H_SAPIENS = {
    "F": {"TTT": 0.45, "TTC": 0.55},
    "L": {"TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40},
    "I": {"ATT": 0.36, "ATC": 0.48, "ATA": 0.16},
    "M": {"ATG": 1.00},
    "V": {"GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47},
    "S": {"TCT": 0.18, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.24},
    "P": {"CCT": 0.28, "CCC": 0.33, "CCA": 0.28, "CCG": 0.11},
    "T": {"ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12},
    "A": {"GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11},
    "Y": {"TAT": 0.43, "TAC": 0.57},
    "H": {"CAT": 0.41, "CAC": 0.59},
    "Q": {"CAA": 0.25, "CAG": 0.75},
    "N": {"AAT": 0.46, "AAC": 0.54},
    "K": {"AAA": 0.42, "AAG": 0.58},
    "D": {"GAT": 0.46, "GAC": 0.54},
    "E": {"GAA": 0.42, "GAG": 0.58},
    "C": {"TGT": 0.45, "TGC": 0.55},
    "W": {"TGG": 1.00},
    "R": {"CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.21},
    "G": {"GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25},
    "*": {"TAA": 0.28, "TAG": 0.20, "TGA": 0.52},
}

# M. musculus is very close to human; bundled as approximate, flagged in note.
M_MUSCULUS = {
    "F": {"TTT": 0.43, "TTC": 0.57},
    "L": {"TTA": 0.06, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.08, "CTG": 0.39},
    "I": {"ATT": 0.34, "ATC": 0.50, "ATA": 0.16},
    "M": {"ATG": 1.00},
    "V": {"GTT": 0.17, "GTC": 0.25, "GTA": 0.12, "GTG": 0.46},
    "S": {"TCT": 0.19, "TCC": 0.22, "TCA": 0.14, "TCG": 0.05, "AGT": 0.15, "AGC": 0.25},
    "P": {"CCT": 0.30, "CCC": 0.31, "CCA": 0.28, "CCG": 0.10},
    "T": {"ACT": 0.25, "ACC": 0.35, "ACA": 0.29, "ACG": 0.11},
    "A": {"GCT": 0.29, "GCC": 0.38, "GCA": 0.23, "GCG": 0.10},
    "Y": {"TAT": 0.42, "TAC": 0.58},
    "H": {"CAT": 0.40, "CAC": 0.60},
    "Q": {"CAA": 0.25, "CAG": 0.75},
    "N": {"AAT": 0.43, "AAC": 0.57},
    "K": {"AAA": 0.39, "AAG": 0.61},
    "D": {"GAT": 0.45, "GAC": 0.55},
    "E": {"GAA": 0.41, "GAG": 0.59},
    "C": {"TGT": 0.48, "TGC": 0.52},
    "W": {"TGG": 1.00},
    "R": {"CGT": 0.09, "CGC": 0.18, "CGA": 0.12, "CGG": 0.19, "AGA": 0.21, "AGG": 0.21},
    "G": {"GGT": 0.18, "GGC": 0.33, "GGA": 0.26, "GGG": 0.23},
    "*": {"TAA": 0.28, "TAG": 0.21, "TGA": 0.51},
}

E_COLI = {
    "F": {"TTT": 0.58, "TTC": 0.42},
    "L": {"TTA": 0.14, "TTG": 0.13, "CTT": 0.12, "CTC": 0.10, "CTA": 0.04, "CTG": 0.47},
    "I": {"ATT": 0.49, "ATC": 0.39, "ATA": 0.12},
    "M": {"ATG": 1.00},
    "V": {"GTT": 0.28, "GTC": 0.20, "GTA": 0.17, "GTG": 0.35},
    "S": {"TCT": 0.15, "TCC": 0.15, "TCA": 0.12, "TCG": 0.15, "AGT": 0.15, "AGC": 0.28},
    "P": {"CCT": 0.18, "CCC": 0.13, "CCA": 0.20, "CCG": 0.49},
    "T": {"ACT": 0.19, "ACC": 0.40, "ACA": 0.16, "ACG": 0.25},
    "A": {"GCT": 0.18, "GCC": 0.26, "GCA": 0.23, "GCG": 0.33},
    "Y": {"TAT": 0.59, "TAC": 0.41},
    "H": {"CAT": 0.57, "CAC": 0.43},
    "Q": {"CAA": 0.34, "CAG": 0.66},
    "N": {"AAT": 0.49, "AAC": 0.51},
    "K": {"AAA": 0.74, "AAG": 0.26},
    "D": {"GAT": 0.63, "GAC": 0.37},
    "E": {"GAA": 0.68, "GAG": 0.32},
    "C": {"TGT": 0.46, "TGC": 0.54},
    "W": {"TGG": 1.00},
    "R": {"CGT": 0.36, "CGC": 0.36, "CGA": 0.07, "CGG": 0.11, "AGA": 0.07, "AGG": 0.03},
    "G": {"GGT": 0.35, "GGC": 0.37, "GGA": 0.13, "GGG": 0.15},
    "*": {"TAA": 0.61, "TAG": 0.09, "TGA": 0.30},
}

TABLES = {
    "h_sapiens": {
        "name": "Homo sapiens",
        "fractions_by_aa": H_SAPIENS,
        "source": "Approximate relative synonymous codon usage, consistent with the "
        "Kazusa Codon Usage Database (https://www.kazusa.or.jp/codon/) "
        "H. sapiens [gbpri] entry.",
        "version": "bundled-snapshot-2024.1",
        "note": "BUNDLED APPROXIMATION. Verify/refresh against the cited source "
        "before production use; values are rounded per-family fractions.",
    },
    "m_musculus": {
        "name": "Mus musculus",
        "fractions_by_aa": M_MUSCULUS,
        "source": "Approximate relative synonymous codon usage, consistent with the "
        "Kazusa Codon Usage Database M. musculus entry.",
        "version": "bundled-snapshot-2024.1",
        "note": "BUNDLED APPROXIMATION; mammalian usage is close to human. "
        "Verify/refresh against the cited source before production use.",
    },
    "e_coli": {
        "name": "Escherichia coli K-12",
        "fractions_by_aa": E_COLI,
        "source": "Approximate relative synonymous codon usage, consistent with the "
        "Kazusa Codon Usage Database E. coli K-12 entry.",
        "version": "bundled-snapshot-2024.1",
        "note": "BUNDLED APPROXIMATION. Verify/refresh against the cited source "
        "before production use.",
    },
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for key, tbl in TABLES.items():
        # sanity: each family should sum to ~1.0
        for aa, fam in tbl["fractions_by_aa"].items():
            total = sum(fam.values())
            if abs(total - 1.0) > 0.03:
                raise SystemExit(f"{key} {aa} sums to {total:.3f}, expected ~1.0")
        path = OUT / f"{key}.json"
        path.write_text(json.dumps({"key": key, **tbl}, indent=2))
        print("wrote", path)


if __name__ == "__main__":
    main()
