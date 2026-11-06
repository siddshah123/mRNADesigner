"""Transcript element library + transcript assembly.

Elements are versioned and source-backed; each carries an explicit uncertainty
statement. Synthetic placeholder UTRs are labeled as such. User-supplied UTRs
and poly(A)/cap choices are honored. Cap chemistry and modified nucleosides are
metadata only and are never written as letters into the sequence.
"""
from __future__ import annotations

import json
import pathlib
from functools import lru_cache
from typing import Optional

from .schemas import (
    CapSpec,
    ModNucleoside,
    PolyASpec,
    Provenance,
    Transcript,
    TranscriptElement,
)

_DATA = pathlib.Path(__file__).parent / "data" / "references"


@lru_cache(maxsize=1)
def _lib() -> dict:
    return json.loads((_DATA / "transcript_elements.json").read_text())


def library_version() -> str:
    return _lib()["library_version"]


def _mk_element(raw: dict) -> TranscriptElement:
    return TranscriptElement(
        element_id=raw["element_id"],
        role=raw["role"],
        name=raw["name"],
        sequence=raw["sequence"],
        version=raw["version"],
        source_id=raw["element_id"],
        relevance=raw["relevance"],
        uncertainty=raw["uncertainty"],
        provenance=Provenance(raw["provenance"]),
    )


def list_elements(role: str) -> list[TranscriptElement]:
    key = {"5utr": "five_utr", "3utr": "three_utr", "kozak": "kozak"}[role]
    return [_mk_element(r) for r in _lib()[key]]


def get_element(role: str, element_id: str) -> Optional[TranscriptElement]:
    for el in list_elements(role):
        if el.element_id == element_id:
            return el
    return None


def default_element(role: str) -> TranscriptElement:
    return list_elements(role)[0]


def list_caps() -> list[dict]:
    return _lib()["cap"]


def list_mods() -> list[dict]:
    return _lib()["mod_nucleosides"]


def polyA_defaults() -> dict:
    return _lib()["polyA_defaults"]


def user_element(role: str, sequence: str) -> TranscriptElement:
    base = get_element(role, f"user_{role}")
    seq = sequence.strip().upper().replace("T", "U") if role != "kozak" else sequence.strip().upper()
    return TranscriptElement(
        element_id=f"user_{role}",
        role=role,  # type: ignore[arg-type]
        name=f"User-supplied {role}",
        sequence=seq,
        version="user",
        source_id=f"user_{role}",
        relevance=(base.relevance if base else "User-provided element."),
        uncertainty="Provenance and validation are the user's responsibility.",
        provenance=Provenance.USER,
    )


def build_transcript(
    intended_protein: str,
    optimized_cds_rna: str,
    *,
    five_utr: TranscriptElement | None = None,
    kozak: TranscriptElement | None = None,
    three_utr: TranscriptElement | None = None,
    polyA: PolyASpec | None = None,
    cap: CapSpec | None = None,
    mods: list[ModNucleoside] | None = None,
    use_defaults: bool = True,
) -> Transcript:
    if use_defaults:
        five_utr = five_utr or default_element("5utr")
        kozak = kozak or default_element("kozak")
        three_utr = three_utr or default_element("3utr")
        if polyA is None:
            d = polyA_defaults()
            polyA = PolyASpec(
                mode=d["mode"], length=d["length"],
                provenance=Provenance.ASSUMPTION, source_id="polyA_defaults",
                note=d["note"],
            )
        if cap is None:
            c = list_caps()[0]
            cap = CapSpec(name=c["name"], chemistry=c["chemistry"],
                          provenance=Provenance.ASSUMPTION, note=c["note"])
    return Transcript(
        intended_protein=intended_protein,
        optimized_cds=optimized_cds_rna,
        five_utr=five_utr,
        kozak=kozak,
        three_utr=three_utr,
        polyA=polyA or PolyASpec(),
        cap=cap or CapSpec(),
        mod_nucleosides=mods or [],
    )
