"use strict";
const $ = (id) => document.getElementById(id);
const COLORS = ["#4aa8ff", "#2ea043", "#d29922", "#a371f7", "#f85149", "#56d4dd", "#e3b341"];
let META = null;
let LAST_ARCH = null; // editable architecture state

async function api(path, body) {
  const r = await fetch(path, {
    method: "POST", headers: { "content-type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!r.ok) {
    const t = await r.text();
    throw new Error(`${r.status}: ${t}`);
  }
  return r.json();
}

function buildRequest() {
  const ctx = {
    protein_species: $("protein_species").value || null,
    expression_host: $("host").value || null,
    cell_type: $("cell_type").value || null,
    localization: $("localization").value || null,
    delivery: $("delivery").value || null,
    in_vivo: $("in_vivo").value === "" ? null : $("in_vivo").value === "true",
    max_total_nt: null,
  };
  const preferred = {};
  if ($("five_utr").value.trim()) preferred.five_utr = $("five_utr").value.trim();
  if ($("three_utr").value.trim()) preferred.three_utr = $("three_utr").value.trim();
  const req = {
    text: $("text").value,
    protein_identifier: $("identifier").value || null,
    pasted_protein_sequence: $("pasted").value || null,
    user_nucleotide_sequence: null,
    preferred_elements: preferred,
    context: ctx,
    architecture_override: LAST_ARCH, // if edited, send it back verbatim
  };
  return req;
}

function buildObjectives() {
  const motifs = $("motifs").value.split(",").map((s) => s.trim()).filter(Boolean);
  return {
    host: $("host").value || "h_sapiens",
    target_gc_min: parseFloat($("gc_min").value),
    target_gc_max: parseFloat($("gc_max").value),
    gc_window: 40,
    max_homopolymer: parseInt($("max_hp").value, 10),
    forbidden_motifs: motifs,
    avoid_restriction_sites: [],
    weight_codon: 1, weight_gc: 1, weight_homopolymer: 1, weight_motif: 1,
    n_candidates: parseInt($("ncand").value, 10),
    seed: parseInt($("seed").value, 10),
  };
}

function body(skip) {
  return {
    request: buildRequest(),
    objectives: buildObjectives(),
    allow_network: $("allow_network").value === "true",
    size_profile_id: $("size_profile").value || null,
    skip_optimization: !!skip,
  };
}

async function loadMeta() {
  META = await (await fetch("/api/meta")).json();
  $("versions").textContent =
    `codon:${META.codon_tables.length} tables · elements:${META.elements_library_version} · lit:${META.literature_corpus_version}`;
  const host = $("host");
  META.codon_tables.forEach((t) => {
    const o = document.createElement("option");
    o.value = t.key; o.textContent = `${t.name} (${t.key})`;
    host.appendChild(o);
  });
  const sp = $("size_profile");
  META.size_profiles.forEach((p) => {
    const o = document.createElement("option");
    o.value = p.profile_id; o.textContent = p.name;
    sp.appendChild(o);
  });
}

// ---------- rendering ----------
function show(id) { $(id).style.display = "block"; }

function renderStatus(design) {
  show("status-line");
  const b = $("status-badges");
  b.innerHTML = "";
  const badge = (cls, txt) => {
    const s = document.createElement("span");
    s.className = `badge ${cls}`; s.textContent = txt; s.style.marginRight = "6px";
    return s;
  };
  const st = design.status;
  b.appendChild(badge(st === "computationally_validated" ? "validated" : (st === "blocked" ? "blocked" : "provisional"),
    st.replace(/_/g, " ")));
  if (design.is_demo_fixture) b.appendChild(badge("demo", "demo fixture"));
  $("status-note").textContent = design.is_demo_fixture
    ? "This design uses BUNDLED SYNTHETIC demo data, not live research."
    : "";
}

function renderClar(clar) {
  if (!clar.questions.length) { $("clar-panel").style.display = "none"; return; }
  show("clar-panel");
  const el = $("clar-body"); el.innerHTML = "";
  clar.questions.forEach((q) => {
    const d = document.createElement("div");
    d.className = "comp";
    d.innerHTML = `<div><b>${q.question}</b></div><div class="hint">Why it matters: ${q.impact}</div>`;
    if (q.options.length) {
      const sel = document.createElement("select");
      sel.innerHTML = `<option value="">(choose)</option>` +
        q.options.map((o) => `<option value="${o.value}">${o.label}</option>`).join("");
      sel.dataset.key = q.key;
      sel.addEventListener("change", () => applyClarification(q.key, sel.value));
      d.appendChild(sel);
    }
    el.appendChild(d);
  });
}

function applyClarification(key, value) {
  if (!value) return;
  if (key === "species_and_host") {
    const [sp, host] = value.split("/");
    $("protein_species").value = sp;
    $("host").value = host;
  } else if (key === "localization") {
    $("localization").value = value;
  } else if (key === "fusion_orientation" && LAST_ARCH) {
    if (value === "reversed") {
      const prot = LAST_ARCH.fusion_order.filter((i) => ["protein", "domain"].includes(LAST_ARCH.components[i].kind));
      if (prot.length >= 2) {
        const a = LAST_ARCH.fusion_order.indexOf(prot[0]);
        const b = LAST_ARCH.fusion_order.indexOf(prot[1]);
        [LAST_ARCH.fusion_order[a], LAST_ARCH.fusion_order[b]] = [LAST_ARCH.fusion_order[b], LAST_ARCH.fusion_order[a]];
        renderArch(LAST_ARCH);
      }
    }
  }
}

function renderArch(arch) {
  LAST_ARCH = arch;
  show("arch-panel");
  const el = $("arch-body"); el.innerHTML = "";
  const title = document.createElement("div");
  title.innerHTML = `<label>Title</label>`;
  const ti = document.createElement("input"); ti.value = arch.title;
  ti.addEventListener("change", () => { arch.title = ti.value; });
  title.appendChild(ti); el.appendChild(title);

  const order = arch.fusion_order.length ? arch.fusion_order : arch.components.map((_, i) => i);
  el.appendChild(Object.assign(document.createElement("div"),
    { className: "hint", textContent: "Fusion order (N→C). Edit sequences/placement below, then re-run Design." }));

  order.forEach((idx, rank) => {
    const c = arch.components[idx];
    const d = document.createElement("div");
    d.className = "comp" + (c.resolved || !["protein", "domain"].includes(c.kind) ? "" : " unresolved");
    d.style.borderLeft = `4px solid ${COLORS[rank % COLORS.length]}`;
    d.innerHTML =
      `<div class="k">${rank + 1}. ${c.kind}${c.placement ? " · " + c.placement : ""} · provenance=${c.provenance}${c.resolved ? "" : " · <span class='tag-q'>UNRESOLVED</span>"}</div>
       <div><b>${c.name}</b> ${c.accession ? `<span class="pill">${c.accession}${c.accession_version ? " v" + c.accession_version : ""}</span>` : ""}</div>`;
    const seq = document.createElement("textarea");
    seq.rows = c.sequence && c.sequence.length > 60 ? 2 : 1;
    seq.value = c.sequence || "";
    seq.placeholder = "(unresolved — paste amino acids here to resolve)";
    seq.style.marginTop = "5px";
    seq.addEventListener("change", () => {
      c.sequence = seq.value.replace(/[^A-Za-z]/g, "").toUpperCase();
      c.resolved = c.sequence.length > 0;
      if (c.resolved) c.provenance = "user";
    });
    d.appendChild(seq);
    if (c.note) d.appendChild(Object.assign(document.createElement("div"), { className: "hint", textContent: c.note }));
    el.appendChild(d);
  });

  if (arch.assumptions && arch.assumptions.length) {
    const a = document.createElement("div"); a.className = "hint";
    a.innerHTML = "<b>Assumptions:</b> " + arch.assumptions.map((x) => `${x.field}=${x.value} (${x.rationale})`).join("; ");
    el.appendChild(a);
  }
  if (arch.unresolved && arch.unresolved.length) {
    const u = document.createElement("div"); u.className = "err-box";
    u.innerHTML = "<b>Unresolved dependencies:</b><br>" + arch.unresolved.join("<br>");
    el.appendChild(u);
  }
  const reb = document.createElement("button");
  reb.textContent = "Re-run Design with edits"; reb.style.marginTop = "10px";
  reb.addEventListener("click", runDesign);
  el.appendChild(reb);
}

function renderEvidence(design) {
  show("evidence-panel");
  const el = $("evidence-body"); el.innerHTML = "";
  if (!design.evidence.length) el.innerHTML = "<div class='hint'>No evidence records matched.</div>";
  design.evidence.forEach((e) => {
    const d = document.createElement("div"); d.className = "evidence";
    d.innerHTML =
      `<div><span class="pill">${e.match_quality}</span> <b>${e.title}</b></div>
       <div class="hint">${e.doi ? "DOI " + e.doi + " · " : ""}${e.pmid ? "PMID " + e.pmid + " · " : ""}${e.accession ? "acc " + e.accession + " · " : ""}retrieved ${e.retrieved} · provenance=${e.provenance}</div>
       <div>${e.finding}</div>
       <div class="hint">Experimental context: ${e.experimental_context || "n/a"} · exact sequence available: <b>${e.exact_sequence_available}</b> · supports: ${e.supports_decision || "—"}</div>
       ${e.url ? `<div><a href="${e.url}" target="_blank" rel="noopener">${e.url}</a></div>` : ""}`;
    el.appendChild(d);
  });
  const rs = document.createElement("details");
  rs.innerHTML = "<summary>Retrieval states (what was contacted)</summary>" +
    design.retrieval_states.map((r) =>
      `<div class="hint">• ${r.integration}: attempted=${r.attempted} live=${r.live} ok=${r.ok} — ${r.detail}</div>`).join("");
  el.appendChild(rs);
}

function renderMapAndSeq(design) {
  if (!design.transcript) return;
  const t = design.transcript;
  show("map-panel"); show("seq-panel");
  // transcript map
  const segs = [];
  if (t.five_utr) segs.push(["5'UTR", t.five_utr.sequence.length, "#56d4dd"]);
  if (t.kozak) segs.push(["Kozak", t.kozak.sequence.length, "#e3b341"]);
  segs.push(["CDS", t.optimized_cds.length, "#4aa8ff"]);
  if (t.three_utr) segs.push(["3'UTR", t.three_utr.sequence.length, "#a371f7"]);
  if (t.polyA.mode === "encoded") segs.push([`polyA(${t.polyA.length})`, t.polyA.length, "#2ea043"]);
  const total = segs.reduce((s, x) => s + x[1], 0);
  const map = $("map-body"); map.innerHTML = "";
  const bar = document.createElement("div"); bar.className = "map";
  segs.forEach(([name, len, col]) => {
    const s = document.createElement("div");
    s.className = "seg2"; s.style.background = col; s.style.width = (100 * len / total) + "%";
    s.title = `${name}: ${len} nt`; s.textContent = len > total * 0.06 ? name : "";
    bar.appendChild(s);
  });
  map.appendChild(bar);
  const meta = document.createElement("div"); meta.className = "hint";
  meta.innerHTML = `cap=<b>${t.cap.name}</b> (${t.cap.chemistry}) — metadata only, not in sequence · ` +
    `poly(A) mode=<b>${t.polyA.mode}</b> len=${t.polyA.length} · no promoter in mature RNA` +
    (t.mod_nucleosides.length ? ` · mods: ${t.mod_nucleosides.map((m) => m.name).join(", ")} (metadata)` : "");
  map.appendChild(meta);
  if (t.five_utr) map.appendChild(utrNote(t.five_utr));
  if (t.three_utr) map.appendChild(utrNote(t.three_utr));

  // sequence viewer: CDS colored by component (aa->nt), UTR flanks neutral
  const sv = $("seq-body"); sv.innerHTML = "";
  const asm = design._assembly;
  const spans = (asm && asm.spans) ? asm.spans.filter((s) => s.end_aa >= s.start_aa) : [];
  const legend = document.createElement("div"); legend.className = "legend";
  spans.forEach((s, i) => {
    const sp = document.createElement("span");
    sp.innerHTML = `<span class="sw" style="background:${COLORS[i % COLORS.length]}"></span>${s.name} (${s.start_aa}–${s.end_aa} aa)`;
    legend.appendChild(sp);
  });
  sv.appendChild(legend);
  const view = document.createElement("div"); view.className = "seqview";
  const cds = t.optimized_cds;
  // build per-nt color from spans
  const colorAt = new Array(cds.length).fill(null);
  spans.forEach((s, i) => {
    const a = (s.start_aa - 1) * 3, b = s.end_aa * 3;
    for (let k = a; k < b && k < cds.length; k++) colorAt[k] = COLORS[i % COLORS.length];
  });
  let html = "<div class='hint'>Optimized CDS (RNA, includes stop), colored by protein component:</div>";
  let cur = null, buf = "";
  const flush = () => { if (buf) html += `<span class="seg" style="background:${cur ? cur + "33" : "transparent"}">${buf}</span>`; buf = ""; };
  for (let k = 0; k < cds.length; k++) {
    if (colorAt[k] !== cur) { flush(); cur = colorAt[k]; }
    buf += cds[k];
  }
  flush();
  view.innerHTML = html;
  sv.appendChild(view);
  const pv = document.createElement("details");
  pv.innerHTML = `<summary>Intended protein (${t.intended_protein.length} aa) & full mature RNA (${mature(t).length} nt)</summary>
    <div class="seqview">${t.intended_protein}</div>
    <div class="hint">Full mature mRNA (5'UTR+Kozak+CDS+3'UTR+encoded polyA; no promoter, cap=metadata):</div>
    <div class="seqview">${mature(t)}</div>`;
  sv.appendChild(pv);
}

function mature(t) {
  let s = "";
  if (t.five_utr) s += t.five_utr.sequence;
  if (t.kozak) s += t.kozak.sequence;
  s += t.optimized_cds;
  if (t.three_utr) s += t.three_utr.sequence;
  if (t.polyA.mode === "encoded") s += "A".repeat(t.polyA.length);
  return s.toUpperCase().replace(/T/g, "U");
}

function utrNote(u) {
  const d = document.createElement("div"); d.className = "hint";
  d.innerHTML = `<b>${u.name}</b> [${u.provenance}] — ${u.relevance} <i>Uncertainty: ${u.uncertainty}</i>`;
  return d;
}

function renderCandidates(opt) {
  if (!opt) { $("cand-panel").style.display = "none"; return; }
  show("cand-panel");
  const el = $("cand-body"); el.innerHTML = "";
  const tbl = document.createElement("table");
  tbl.innerHTML = `<tr><th>candidate</th><th>penalty↓</th><th>CAI</th><th>GC</th><th>max hp</th><th>local GC</th><th>motif hits</th><th></th></tr>` +
    opt.candidates.map((c) => {
      const m = c.metrics;
      const sel = c.candidate_id === opt.selected_id;
      return `<tr class="${sel ? "sel" : ""}"><td>${c.candidate_id}${c.is_baseline ? " <span class='pill'>baseline</span>" : ""}</td>
        <td>${c.objective_score}</td><td>${m.cai}</td><td>${m.gc_content}</td><td>${m.max_homopolymer}</td>
        <td>${m.min_local_gc}–${m.max_local_gc}</td><td>${m.forbidden_motif_hits}</td>
        <td>${sel ? "✓ selected" : ""}</td></tr>`;
    }).join("");
  el.appendChild(tbl);
  el.appendChild(Object.assign(document.createElement("div"), { className: "hint", textContent: opt.selection_rationale }));
  el.appendChild(Object.assign(document.createElement("div"), { className: "hint",
    textContent: `optimizer=${opt.optimizer_version} · codon_table=${opt.codon_table_version} · seed=${opt.seed} (reproducible)` }));
}

function renderValidation(design) {
  const v = design.validation;
  if (!v) { $("val-panel").style.display = "none"; return; }
  show("val-panel");
  const el = $("val-body"); el.innerHTML = "";
  const group = (title, items, cls) => {
    const h = document.createElement("div");
    h.innerHTML = `<b>${title} (${items.length})</b>`;
    el.appendChild(h);
    items.forEach((i) => {
      const d = document.createElement("div"); d.className = `vitem ${cls}`;
      d.innerHTML = `<span class="c">${i.check}</span> — ${i.message}`;
      el.appendChild(d);
    });
  };
  const errors = v.items.filter((i) => i.severity === "error");
  const warns = v.items.filter((i) => i.severity === "warning");
  const un = v.items.filter((i) => i.severity === "unevaluated");
  const pass = v.items.filter((i) => i.severity === "pass");
  group("Errors", errors, "error");
  group("Warnings", warns, "warning");
  group("Unevaluated", un, "unevaluated");
  const dp = document.createElement("details");
  dp.innerHTML = `<summary>Passed checks (${pass.length})</summary>` +
    pass.map((i) => `<div class="vitem pass"><span class="c">${i.check}</span> — ${i.message}</div>`).join("");
  el.appendChild(dp);

  if (v.size) {
    const s = v.size;
    const sz = document.createElement("div"); sz.className = "comp"; sz.style.marginTop = "8px";
    sz.innerHTML = `<div class="k">size report</div>
      <div class="kv">
        <span class="lbl">protein</span><span>${s.protein_len_aa} aa</span>
        <span class="lbl">CDS</span><span>${s.cds_len_nt} nt (includes stop: ${s.cds_includes_stop})</span>
        <span class="lbl">5′UTR / 3′UTR</span><span>${s.five_utr_len_nt} / ${s.three_utr_len_nt} nt</span>
        <span class="lbl">poly(A)</span><span>${s.polyA_encoded_len_nt} nt (${s.polyA_mode})</span>
        <span class="lbl">TOTAL RNA</span><span><b>${s.total_rna_len_nt} nt</b></span>
        <span class="lbl">software limit</span><span>${s.software_input_limit_nt} nt (not biological)</span>
        <span class="lbl">biological guidance</span><span>${s.biological_guidance} <i>[source=${s.guidance_source_id}, uncertain=${s.guidance_is_uncertain}]</i></span>
      </div>`;
    el.appendChild(sz);
  }
  if (v.translated_protein) {
    const tp = document.createElement("details");
    tp.innerHTML = `<summary>Exact translated protein (${v.translated_protein.length} aa)</summary><div class="seqview">${v.translated_protein}</div>`;
    el.appendChild(tp);
  }
}

function renderDownloads(design, roundtrip) {
  if (!design.transcript) { $("dl-panel").style.display = "none"; return; }
  show("dl-panel");
  const el = $("dl-body"); el.innerHTML = "";
  const formats = [
    ["rna_fasta", "RNA FASTA"], ["cds_fasta", "CDS FASTA"], ["protein_fasta", "Protein FASTA"],
    ["genbank", "GenBank (.gb)"], ["genbank_ivt_template", "IVT DNA template (.gb)"],
    ["manifest", "JSON manifest"], ["report", "Report (.txt)"],
  ];
  formats.forEach(([fmt, label]) => {
    const b = document.createElement("button"); b.className = "secondary"; b.textContent = label;
    b.addEventListener("click", () => downloadExport(fmt, label));
    el.appendChild(b);
  });
  if (roundtrip) {
    $("roundtrip").innerHTML =
      `Export round-trip (BioPython re-parse): GenBank=<b>${roundtrip.genbank_roundtrip}</b>, FASTA=<b>${roundtrip.fasta_roundtrip}</b> — ${roundtrip.detail}`;
  }
}

async function downloadExport(fmt, label) {
  const r = await fetch(`/api/export/${fmt}`, {
    method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(body(false)),
  });
  if (!r.ok) { alert(`Export failed: ${await r.text()}`); return; }
  const text = await r.text();
  const ext = fmt === "manifest" ? "json" : (fmt.includes("genbank") ? "gb" : (fmt === "report" ? "txt" : "fasta"));
  const blob = new Blob([text], { type: "text/plain" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = `${fmt}.${ext}`;
  a.click(); URL.revokeObjectURL(a.href);
}

// ---------- actions ----------
async function runInterpret() {
  LAST_ARCH = null;
  try {
    const res = await api("/api/interpret", body(true));
    renderClar(res.clarification);
    renderArch(res.architecture);
    ["status-line", "evidence-panel", "map-panel", "seq-panel", "cand-panel", "val-panel", "dl-panel"]
      .forEach((id) => $(id).style.display = "none");
  } catch (e) { alert(e.message); }
}

function hasInput() {
  return ($("text").value.trim() || $("identifier").value.trim() || $("pasted").value.trim());
}

async function runDesign() {
  if (!hasInput()) { alert("Enter at least a description in box 1 (the only required field)."); return; }
  const btn = $("btn-design"); const old = btn.innerHTML;
  btn.innerHTML = '<span class="spin"></span> designing…'; btn.disabled = true;
  try {
    const res = await api("/api/design", body(false));
    const design = res.design;
    design._assembly = await fetchAssembly(); // for sequence highlighting
    renderStatus(design);
    renderClar(res.clarification);
    renderArch(design.architecture);
    renderEvidence(design);
    renderMapAndSeq(design);
    renderCandidates(design.optimization);
    renderValidation(design);
    renderDownloads(design, res.roundtrip);
  } catch (e) {
    alert(e.message);
  } finally {
    btn.innerHTML = old; btn.disabled = false;
  }
}

// assembly spans are derivable from architecture; ask interpret endpoint isn't enough,
// so we reconstruct them client-side from component lengths in fusion order.
async function fetchAssembly() {
  if (!LAST_ARCH) return null;
  const order = LAST_ARCH.fusion_order.length ? LAST_ARCH.fusion_order : LAST_ARCH.components.map((_, i) => i);
  const spans = []; let cursor = 0;
  order.forEach((idx) => {
    const c = LAST_ARCH.components[idx];
    if (["protein", "domain"].includes(c.kind) && (!c.resolved || !c.sequence)) return;
    const len = (c.sequence || "").length;
    if (len === 0) return;
    spans.push({ name: c.name, kind: c.kind, start_aa: cursor + 1, end_aa: cursor + len });
    cursor += len;
  });
  return { spans };
}

const EXAMPLES = {
  tagged: { identifier: "demo_reporter", text: "demo reporter with a C-terminal His tag" },
  fusion: { identifier: "demo_reporter fused to demo_secreted_precursor", text: "two-domain fusion with a linker" },
  precursor: { identifier: "demo_secreted_precursor", text: "secreted demo precursor", localization: "secreted" },
  ambiguous: { identifier: "interleukin-2", text: "a cytokine" },
  toolong: { identifier: "demo_reporter", text: "length test", _long5: true },
};
function applyExample(key) {
  const ex = EXAMPLES[key];
  $("identifier").value = ex.identifier || "";
  $("text").value = ex.text || "";
  $("pasted").value = "";
  $("localization").value = ex.localization || "";
  $("size_profile").value = key === "toolong" ? "lnp_ivt_mrna" : "";
  $("five_utr").value = ex._long5 ? "A".repeat(6000) : "";
  $("three_utr").value = "";
}

window.addEventListener("DOMContentLoaded", async () => {
  await loadMeta();
  $("btn-interpret").addEventListener("click", runInterpret);
  $("btn-design").addEventListener("click", runDesign);
  document.querySelectorAll("[data-ex]").forEach((a) =>
    a.addEventListener("click", (e) => { e.preventDefault(); applyExample(a.dataset.ex); }));
});
