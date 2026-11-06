"""Reference-sequence retrieval.

Two paths, both honest about provenance:
  * Bundled library: standard peptide tags + linkers (verifiable common-knowledge
    design elements) and clearly-labeled SYNTHETIC demo proteins.
  * Live UniProt REST fetch for real accessions. Records a RetrievalState with
    (attempted, live, ok, detail). On failure it does NOT fabricate a sequence;
    callers mark the affected component unresolved/provisional.
"""
from __future__ import annotations

import json
import pathlib
import re
from datetime import date
from functools import lru_cache
from typing import Optional

from pydantic import BaseModel

from .schemas import Provenance, RetrievalState

_DATA = pathlib.Path(__file__).parent / "data" / "references"


class ReferenceSequence(BaseModel):
    name: str
    sequence: str  # amino acids
    species: Optional[str] = None
    accession: Optional[str] = None
    accession_version: Optional[str] = None
    is_precursor: bool = False
    signal_peptide: Optional[str] = None
    mature_sequence: Optional[str] = None
    mature_note: str = ""
    source: str = ""
    provenance: Provenance = Provenance.REFERENCE
    is_fixture: bool = False
    retrieval: Optional[RetrievalState] = None


@lru_cache(maxsize=1)
def _elements() -> dict:
    return json.loads((_DATA / "elements.json").read_text())


def library_version() -> str:
    return _elements()["library_version"]


def list_tags() -> list[dict]:
    return _elements()["tags"]


def list_linkers() -> list[dict]:
    return _elements()["linkers"]


def list_demo_proteins() -> list[dict]:
    return _elements()["demo_proteins"]


def get_tag(name: str) -> Optional[dict]:
    for t in list_tags():
        if t["name"].lower() == name.lower():
            return t
    return None


def get_linker(name: str) -> Optional[dict]:
    for lk in list_linkers():
        if lk["name"].lower() == name.lower():
            return lk
    return None


def get_demo_protein(key: str) -> Optional[dict]:
    for p in list_demo_proteins():
        if p["key"].lower() == key.lower() or p["accession"].lower() == key.lower():
            return p
    return None


_ACCESSION_RE = re.compile(r"^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")


def looks_like_uniprot_accession(s: str) -> bool:
    return bool(_ACCESSION_RE.match(s.strip().upper()))


def fetch_uniprot(accession: str, allow_network: bool = True, timeout: float = 20.0) -> ReferenceSequence:
    """Attempt a live UniProt fetch. Always returns a ReferenceSequence whose
    `.retrieval` records exactly what happened. On any failure, `.sequence` is
    empty and provenance stays unresolved (never fabricated)."""
    acc = accession.strip().upper()
    state = RetrievalState(integration="uniprot", attempted=True, live=False, ok=False)

    if not allow_network:
        state.detail = "network disabled by caller (allow_network=False)"
        return ReferenceSequence(name=acc, sequence="", accession=acc,
                                 provenance=Provenance.ASSUMPTION, retrieval=state)
    try:
        import httpx

        state.live = True
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        r = httpx.get(url, timeout=timeout, follow_redirects=True)
        if r.status_code != 200 or not r.text.startswith(">"):
            state.detail = f"HTTP {r.status_code} from UniProt for {acc}"
            return ReferenceSequence(name=acc, sequence="", accession=acc,
                                     provenance=Provenance.ASSUMPTION, retrieval=state)
        lines = r.text.splitlines()
        header = lines[0]
        seq = "".join(lines[1:]).replace(" ", "").upper()
        # parse species + version from header if present
        species = None
        m = re.search(r"OS=([^=]+?)(?: [A-Z]{2}=|$)", header)
        if m:
            species = m.group(1).strip()
        ver = None
        m2 = re.search(r"SV=(\d+)", header)
        if m2:
            ver = m2.group(1)
        state.ok = True
        state.detail = f"fetched {len(seq)} aa from UniProt {acc}"
        return ReferenceSequence(
            name=acc, sequence=seq, species=species, accession=acc,
            accession_version=ver, source=f"UniProt {acc} ({url})",
            provenance=Provenance.REFERENCE, retrieval=state,
        )
    except Exception as e:  # network blocked, timeout, etc.
        state.detail = f"{type(e).__name__}: {str(e)[:160]}"
        return ReferenceSequence(name=acc, sequence="", accession=acc,
                                 provenance=Provenance.ASSUMPTION, retrieval=state)


def search_uniprot(
    name: str,
    species: str | None = None,
    allow_network: bool = True,
    timeout: float = 25.0,
) -> ReferenceSequence:
    """Resolve a protein NAME to a real sequence via the UniProt search API
    (deterministic 'research' — no LLM). Prefers reviewed (Swiss-Prot) entries,
    optionally filtered by organism. Returns the top hit. On any failure the
    sequence is empty and the retrieval state records why (never fabricated).

    Requires outbound access to rest.uniprot.org (works on a normal machine;
    blocked inside network-restricted sandboxes)."""
    q = name.strip()
    state = RetrievalState(integration="uniprot_search", attempted=True, live=False, ok=False)
    if not q:
        state.detail = "empty query"
        return ReferenceSequence(name=name, sequence="", provenance=Provenance.ASSUMPTION, retrieval=state)
    if not allow_network:
        state.detail = ("network disabled (allow_network=False); enable 'Network retrieval' "
                        "to let UniProt resolve this protein name.")
        return ReferenceSequence(name=name, sequence="", provenance=Provenance.ASSUMPTION, retrieval=state)
    try:
        import httpx

        state.live = True
        # Build a UniProt query string. Reviewed entries first; optional organism.
        query = f'("{q}") AND reviewed:true'
        if species:
            query += f' AND organism_name:"{species}"'
        params = {
            "query": query,
            "format": "fasta",
            "size": "1",
            "sort": "annotation_score desc",
        }
        url = "https://rest.uniprot.org/uniprotkb/search"
        r = httpx.get(url, params=params, timeout=timeout, follow_redirects=True)
        if r.status_code != 200 or not r.text.startswith(">"):
            state.detail = f"HTTP {r.status_code} or no reviewed hit for {q!r}"
            return ReferenceSequence(name=name, sequence="", provenance=Provenance.ASSUMPTION, retrieval=state)
        lines = r.text.splitlines()
        header = lines[0]
        seq = "".join(l for l in lines[1:] if not l.startswith(">")).replace(" ", "").upper()
        # header: >sp|ACCESSION|NAME ... OS=Species ... SV=n
        acc = None
        hm = re.match(r">\w+\|([^|]+)\|", header)
        if hm:
            acc = hm.group(1)
        sp = None
        m = re.search(r"OS=([^=]+?)(?: [A-Z]{2}=|$)", header)
        if m:
            sp = m.group(1).strip()
        ver = None
        m2 = re.search(r"SV=(\d+)", header)
        if m2:
            ver = m2.group(1)
        state.ok = True
        state.detail = f"UniProt search {q!r} -> {acc} ({len(seq)} aa)"
        return ReferenceSequence(
            name=f"{q} (UniProt {acc})", sequence=seq, species=sp, accession=acc,
            accession_version=ver,
            source=f"UniProt search for {q!r}; top reviewed hit {acc} "
            f"(https://www.uniprot.org/uniprotkb/{acc}).",
            provenance=Provenance.REFERENCE, retrieval=state,
        )
    except Exception as e:
        state.detail = f"{type(e).__name__}: {str(e)[:160]}"
        return ReferenceSequence(name=name, sequence="", provenance=Provenance.ASSUMPTION, retrieval=state)


def demo_as_reference(key: str) -> Optional[ReferenceSequence]:
    p = get_demo_protein(key)
    if not p:
        return None
    return ReferenceSequence(
        name=p["display_name"],
        sequence=p["sequence"],
        species=p.get("species"),
        accession=p.get("accession"),
        is_precursor=p.get("is_precursor", False),
        signal_peptide=p.get("signal_peptide"),
        mature_sequence=p.get("mature_sequence"),
        mature_note=p.get("mature_note", ""),
        source=p.get("source", ""),
        provenance=Provenance.FIXTURE,
        is_fixture=True,
        retrieval=RetrievalState(integration="bundled_fixture", attempted=True,
                                 live=False, ok=True, detail="bundled synthetic demo"),
    )
