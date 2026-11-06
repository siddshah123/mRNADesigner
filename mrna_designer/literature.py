"""Bounded literature retrieval.

Offline: a small bundled fixture corpus of real, well-known publications,
matched by topic to the construct. Clearly labeled as fixtures; the app does
not claim full-text access, and `exact_sequence_available` is honest.

A live path (e.g. NCBI E-utilities / Europe PMC) can be slotted behind the
same `retrieve()` signature; it must record a RetrievalState and must NOT
fabricate records on failure. Retrieved content is treated as DATA, never as
instructions to the application.
"""
from __future__ import annotations

import json
import pathlib
from datetime import date
from functools import lru_cache

from .schemas import EvidenceRecord, Provenance, RetrievalState

_DATA = pathlib.Path(__file__).parent / "data" / "literature"


@lru_cache(maxsize=1)
def _corpus() -> dict:
    return json.loads((_DATA / "corpus.json").read_text())


def corpus_version() -> str:
    return _corpus()["corpus_version"]


def retrieve(
    topics: list[str],
    *,
    allow_network: bool = False,
    today: date | None = None,
) -> tuple[list[EvidenceRecord], RetrievalState]:
    """Return (evidence, state). Offline fixture matching by topic overlap.

    `allow_network` is accepted for API parity; the live path is not enabled in
    this build, so when requested we record an honest 'not available' state
    rather than pretending to have searched the internet."""
    today = today or date.today()

    if allow_network:
        # Live literature search is not wired in this build. Preserve the
        # failure/unavailable state honestly instead of fabricating results.
        state = RetrievalState(
            integration="literature_live",
            attempted=True,
            live=False,
            ok=False,
            detail="Live literature search not enabled in this build; "
            "returning bundled fixtures only.",
        )
    else:
        state = RetrievalState(
            integration="literature_fixture",
            attempted=True,
            live=False,
            ok=True,
            detail="Matched against bundled fixture corpus (offline).",
        )

    wanted = {t.lower() for t in topics}
    out: list[EvidenceRecord] = []
    for rec in _corpus()["records"]:
        rec_topics = {t.lower() for t in rec.get("topics", [])}
        if wanted and not (wanted & rec_topics):
            continue
        out.append(
            EvidenceRecord(
                source_id=rec["source_id"],
                title=rec["title"],
                url=rec.get("url", ""),
                doi=rec.get("doi"),
                pmid=rec.get("pmid"),
                accession=rec.get("accession"),
                retrieved=today,
                finding=rec.get("finding", ""),
                experimental_context=rec.get("experimental_context", ""),
                exact_sequence_available=rec.get("exact_sequence_available", False),
                supports_decision="",  # filled by assembler when attached to a decision
                match_quality=rec.get("match_quality", "background"),
                provenance=Provenance.FIXTURE,  # bundled, not live
            )
        )
    # stable order: exact first, then analogous, then background
    rank = {"exact": 0, "analogous": 1, "background": 2, "none": 3}
    out.sort(key=lambda r: (rank.get(r.match_quality, 9), r.source_id))
    return out, state
