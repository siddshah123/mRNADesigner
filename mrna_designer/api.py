"""FastAPI application wiring the pipeline + exports to a browser UI.

Privacy posture:
  * API keys / any external credentials stay server-side (none are required for
    the offline default build).
  * Unpublished user sequences are NOT sent to external services unless the
    request explicitly sets allow_network=True (opt-in), and the response
    discloses every retrieval integration that was contacted (retrieval_states).
"""
from __future__ import annotations

import pathlib

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import codon_tables, export, literature
from . import references as R
from . import transcript as T
from .interpret import interpret
from .pipeline import run_pipeline
from .schemas import ConstructRequest, OptimizationObjectives
from .validate import SIZE_PROFILES

app = FastAPI(title="mRNA Designer", version="0.1.0")

_FRONTEND = pathlib.Path(__file__).parent.parent / "frontend"


class DesignRequestBody(BaseModel):
    request: ConstructRequest
    objectives: OptimizationObjectives | None = None
    allow_network: bool = False
    size_profile_id: str | None = None
    skip_optimization: bool = False


@app.get("/api/meta")
def meta() -> dict:
    return {
        "codon_tables": [
            {"key": k, **{kk: vv for kk, vv in _table_meta(k).items()}}
            for k in codon_tables.available()
        ],
        "tags": R.list_tags(),
        "linkers": R.list_linkers(),
        "demo_proteins": [
            {"key": p["key"], "display_name": p["display_name"],
             "accession": p["accession"], "is_precursor": p.get("is_precursor", False)}
            for p in R.list_demo_proteins()
        ],
        "five_utr": [e.model_dump() for e in T.list_elements("5utr")],
        "three_utr": [e.model_dump() for e in T.list_elements("3utr")],
        "kozak": [e.model_dump() for e in T.list_elements("kozak")],
        "caps": T.list_caps(),
        "mods": T.list_mods(),
        "size_profiles": [p.model_dump() for p in SIZE_PROFILES.values()],
        "elements_library_version": R.library_version(),
        "transcript_library_version": T.library_version(),
        "literature_corpus_version": literature.corpus_version(),
    }


def _table_meta(key: str) -> dict:
    t = codon_tables.load_table(key)
    return {"name": t.name, "version": t.version, "source": t.source, "note": t.note}


@app.post("/api/interpret")
def api_interpret(body: DesignRequestBody) -> dict:
    arch, clar = interpret(body.request, allow_network=body.allow_network)
    return {"architecture": arch.model_dump(), "clarification": clar.model_dump()}


@app.post("/api/design")
def api_design(body: DesignRequestBody) -> dict:
    try:
        design, clar = run_pipeline(
            body.request,
            objectives=body.objectives,
            allow_network=body.allow_network,
            size_profile_id=body.size_profile_id,
            skip_optimization=body.skip_optimization,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    roundtrip = export.roundtrip_check(design) if design.transcript else {
        "genbank_roundtrip": "unevaluated", "fasta_roundtrip": "unevaluated",
        "detail": "no transcript",
    }
    return {
        "design": design.model_dump(),
        "clarification": clar.model_dump(),
        "roundtrip": roundtrip,
    }


# ---- exports: re-run deterministically from the posted request -----------
def _rebuild(body: DesignRequestBody):
    design, _ = run_pipeline(
        body.request, objectives=body.objectives,
        allow_network=body.allow_network, size_profile_id=body.size_profile_id,
    )
    if not design.transcript:
        raise HTTPException(status_code=400, detail="Design is blocked/provisional with no transcript to export.")
    return design


@app.post("/api/export/{fmt}", response_class=PlainTextResponse)
def api_export(fmt: str, body: DesignRequestBody) -> PlainTextResponse:
    design = _rebuild(body)
    try:
        if fmt == "rna_fasta":
            return PlainTextResponse(export.rna_fasta(design))
        if fmt == "cds_fasta":
            return PlainTextResponse(export.cds_fasta(design))
        if fmt == "protein_fasta":
            return PlainTextResponse(export.protein_fasta(design))
        if fmt == "genbank":
            return PlainTextResponse(export.genbank(design))
        if fmt == "genbank_ivt_template":
            return PlainTextResponse(export.genbank(design, as_dna_template=True))
        if fmt == "manifest":
            return PlainTextResponse(export.design_manifest(design), media_type="application/json")
        if fmt == "report":
            return PlainTextResponse(export.report_text(design))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    raise HTTPException(status_code=404, detail=f"unknown export format {fmt!r}")


# ---- static frontend ------------------------------------------------------
@app.get("/", response_class=HTMLResponse)
def index() -> HTMLResponse:
    html = (_FRONTEND / "index.html").read_text()
    return HTMLResponse(html)


if _FRONTEND.exists():
    app.mount("/static", StaticFiles(directory=str(_FRONTEND)), name="static")
