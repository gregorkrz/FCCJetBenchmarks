"""Stage 4 of the resolution pipeline (no ROOT, no fitting — pure packaging).

Reads the consolidated `dashboard_data.json` produced by
`build_dashboard_data.py` and writes a single self-contained HTML file
(Plotly.js loaded from CDN, vanilla JS, no build step) with two tabs:

Explorer tab:
- Multi-select filters for method (jet clustering algo / detector config),
  process, and quantity (jet-part energy resolution, angular resolution,
  eta/cos(theta) scans, or Higgs mass).
- The Higgs mass quantity plots the reco mH histogram (with its Gaussian
  peak fit) per selected method/process, or all four mH definitions
  (reco/gen/GT/reco-GT matched) overlaid when a single combination is
  selected.
- A color picker per selected (method, process) combination, defaulting to
  the colors from `process_config.py` with automatic fallback colors.
- A main Plotly scatter+line plot of resolution vs. energy (or eta/costheta)
  for every selected combination, including the fitted curve.
- Clicking on a data point toggles its underlying per-bin histogram (and the
  low/high/MPV fit values as vertical lines) on/off in a secondary plot, so
  multiple histograms (e.g. from different methods/processes/energy bins)
  can be compared side by side.
- One-click presets (mirroring the comparisons in joint_plots.py and
  presentation.pdf) for common combinations: jet multiplicity (2/4/6 jets),
  clustering algorithm scan (Durham vs. anti-kt radii), detector/matching
  comparison (PF vs. Calo vs. ideal matching), energy recovery on/off,
  B-hadron content scan, and two Higgs-mass peak comparisons (clustering
  algorithm scan, detector comparison) — each with sensible default colors
  that can still be overridden afterwards.
- A "Download raw JSON data" button that re-offers the page's inlined
  DASHBOARD_DATA as a downloadable dashboard_data.json.

Statistics tab: fit coefficients (JER/angular + Higgs mass peak), raw event
counts, and filter pass rates, one table each, built from the same data.

Usage:
    python src/plotting/make_interactive_dashboard.py --inputDir $PATH_TO_HISTOGRAMS/plots/dashboard_data.json
    python src/plotting/make_interactive_dashboard.py --data $PATH_TO_HISTOGRAMS/plots/dashboard_data.json --output dashboard.html
"""
import argparse
import json
import os

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FCC Jet Resolution Dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  body { font-family: sans-serif; margin: 0; padding: 16px; background: #fafafa; }
  h1 { font-size: 20px; margin: 0 0 12px 0; }
  #layout { display: flex; gap: 16px; align-items: flex-start; }
  #controls { flex: 0 0 300px; background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 12px; }
  #plots { flex: 1; min-width: 0; }
  fieldset { border: 1px solid #ccc; border-radius: 4px; margin-bottom: 10px; }
  legend { font-weight: bold; font-size: 12px; color: #444; }
  .scroll-box { max-height: 160px; overflow-y: auto; }
  label.opt { display: block; font-size: 13px; padding: 2px 0; cursor: pointer; }
  #selection-list { font-size: 12px; }
  .sel-row { display: flex; align-items: center; gap: 6px; padding: 2px 0; }
  .sel-row input[type=color] { width: 28px; height: 20px; padding: 0; border: none; }
  .sel-row span { flex: 1; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  #mainPlot { width: 100%; height: 560px; }
  #histPlot { width: 100%; height: 360px; margin-top: 12px; }
  #histHint { font-size: 12px; color: #777; padding: 4px 0; display: flex; align-items: center; gap: 12px; justify-content: space-between; }
  #histSelectionList { font-size: 12px; margin-top: 4px; }
  .hist-row { display: flex; align-items: center; gap: 6px; padding: 2px 0; }
  .hist-row .swatch { width: 12px; height: 12px; border-radius: 2px; flex: 0 0 auto; }
  .hist-row span { flex: 1; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  .hist-row button { border: none; background: #eee; border-radius: 3px; cursor: pointer; font-size: 11px; padding: 1px 6px; }
  #clearHistBtn { border: none; background: #eee; border-radius: 3px; cursor: pointer; font-size: 11px; padding: 2px 8px; }
  #presetList { display: flex; flex-direction: column; gap: 6px; }
  .preset-btn { text-align: left; border: 1px solid #ccc; background: #f5f5f9; border-radius: 4px;
                padding: 6px 8px; cursor: pointer; font-size: 12px; }
  .preset-btn:hover { background: #eaeaf5; }
  .preset-btn .preset-desc { display: block; font-size: 10px; color: #777; font-weight: normal; margin-top: 2px; }
  #header-row { display: flex; align-items: center; justify-content: space-between; gap: 12px; margin-bottom: 8px; }
  #downloadJsonBtn { border: 1px solid #ccc; background: #fff; border-radius: 4px; padding: 6px 10px;
                     cursor: pointer; font-size: 12px; white-space: nowrap; }
  #downloadJsonBtn:hover { background: #eee; }
  #tabBar { display: flex; gap: 4px; margin-bottom: 12px; border-bottom: 1px solid #ccc; }
  .tab-btn { border: 1px solid #ccc; border-bottom: none; background: #eee; border-radius: 6px 6px 0 0;
             padding: 8px 16px; cursor: pointer; font-size: 13px; position: relative; top: 1px; }
  .tab-btn.active { background: #fff; font-weight: bold; border-bottom: 1px solid #fff; }
  .tab-content { display: none; }
  .tab-content.active { display: block; }
  .stats-section { background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 12px; margin-bottom: 16px; }
  .stats-section h2 { font-size: 15px; margin: 0 0 8px 0; }
  .stats-table-wrap { overflow-x: auto; }
  table.stats-table { border-collapse: collapse; font-size: 12px; width: 100%; }
  table.stats-table th, table.stats-table td { border: 1px solid #ddd; padding: 4px 8px; text-align: right; white-space: nowrap; }
  table.stats-table th:first-child, table.stats-table td:first-child { text-align: left; }
  table.stats-table th { background: #f5f5f5; }
  table.stats-table td.best { font-weight: bold; }
</style>
</head>
<body>
<div id="header-row">
  <h1>FCC Jet Resolution Dashboard</h1>
  <button id="downloadJsonBtn">Download raw JSON data</button>
</div>
<div id="tabBar">
  <button class="tab-btn active" data-tab="explorer">Explorer</button>
  <button class="tab-btn" data-tab="statistics">Statistics</button>
</div>

<div id="tab-explorer" class="tab-content active">
<div id="layout">
  <div id="controls">
    <fieldset>
      <legend>Presets</legend>
      <div id="presetList"></div>
    </fieldset>
    <fieldset>
      <legend>Quantity</legend>
      <select id="quantitySelect" style="width: 100%;"></select>
    </fieldset>
    <fieldset>
      <legend>Methods</legend>
      <div class="scroll-box" id="methodOptions"></div>
    </fieldset>
    <fieldset>
      <legend>Processes</legend>
      <div class="scroll-box" id="processOptions"></div>
    </fieldset>
    <fieldset>
      <legend>Selected curves (click swatch to recolor)</legend>
      <div class="scroll-box" id="selection-list"></div>
    </fieldset>
  </div>
  <div id="plots">
    <div id="mainPlot"></div>
    <div id="histHint">
      <span>Click on a data point above to toggle its histogram on/off below. Click multiple points to compare them.</span>
      <label style="white-space:nowrap;"><input type="checkbox" id="normalizeHistCheckbox" checked> normalize</label>
      <button id="clearHistBtn">Clear all</button>
    </div>
    <div id="histPlot"></div>
    <div id="histSelectionList"></div>
  </div>
</div>
</div>

<div id="tab-statistics" class="tab-content">
  <div class="stats-section">
    <h2>Fit coefficients</h2>
    <div class="stats-table-wrap" id="fitCoeffTable"></div>
  </div>
  <div class="stats-section">
    <h2>Event counts (before / after the fully-matched filter)</h2>
    <div class="stats-table-wrap" id="eventCountTable"></div>
  </div>
  <div class="stats-section">
    <h2>Filter pass rates</h2>
    <div class="stats-table-wrap" id="passRateTable"></div>
  </div>
</div>

<script>
const DASHBOARD_DATA = __DASHBOARD_DATA_JSON__;

const QUANTITIES = [
  {key: "energy:_all", label: "Jet energy resolution (all)", kind: "energy", part: "_all"},
  {key: "energy:_charged", label: "Jet energy resolution (charged)", kind: "energy", part: "_charged"},
  {key: "energy:_neutral", label: "Jet energy resolution (neutral)", kind: "energy", part: "_neutral"},
  {key: "energy:_photons", label: "Jet energy resolution (photons)", kind: "energy", part: "_photons"},
  {key: "angle:theta", label: "Angular resolution (theta)", kind: "angle", part: "theta"},
  {key: "angle:phi", label: "Angular resolution (phi)", kind: "angle", part: "phi"},
  {key: "angle:eta", label: "Angular resolution (eta)", kind: "angle", part: "eta"},
  {key: "eta_scan", label: "Energy resolution vs. eta", kind: "eta_scan", part: null},
  {key: "costheta_scan", label: "Energy resolution vs. cos(theta)", kind: "costheta_scan", part: null},
  {key: "mass", label: "Higgs mass (reco / gen / GT / reco-GT matched)", kind: "mass", part: null},
];

const MASS_DEFINITIONS = [
  {key: "reco", label: "reco", color: "#1f77b4"},
  {key: "gen", label: "gen", color: "#ff7f0e"},
  {key: "gt", label: "GT", color: "#2ca02c"},
  {key: "gt_recomatched", label: "reco-GT matched", color: "#d62728"},
];

const AUTO_COLORS = [
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
];

const methodNames = Object.keys(DASHBOARD_DATA.methods).sort();
const processNames = Object.keys(DASHBOARD_DATA.process_meta).sort();

// selection state: Map("method||process" -> color)
const selectedColors = new Map();
let autoColorIdx = 0;

function buildOptionList(container, names, getLabel, onChange) {
  container.innerHTML = "";
  names.forEach(name => {
    const id = "opt-" + name.replace(/[^a-zA-Z0-9]/g, "_") + "-" + Math.random().toString(36).slice(2,7);
    const label = document.createElement("label");
    label.className = "opt";
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.value = name;
    cb.id = id;
    cb.addEventListener("change", onChange);
    label.appendChild(cb);
    label.appendChild(document.createTextNode(" " + getLabel(name)));
    container.appendChild(label);
  });
}

function getChecked(container) {
  return Array.from(container.querySelectorAll("input[type=checkbox]:checked")).map(cb => cb.value);
}

const quantitySelect = document.getElementById("quantitySelect");
QUANTITIES.forEach(q => {
  const opt = document.createElement("option");
  opt.value = q.key;
  opt.textContent = q.label;
  quantitySelect.appendChild(opt);
});

const methodOptions = document.getElementById("methodOptions");
const processOptions = document.getElementById("processOptions");

function getEntry(method, process, quantity) {
  const procData = DASHBOARD_DATA.methods[method].processes[process];
  if (!procData) return null;
  if (quantity.kind === "energy") {
    return (procData.energy || {})[quantity.part] || null;
  }
  if (quantity.kind === "angle") {
    return (procData.angles || {})[quantity.part] || null;
  }
  if (quantity.kind === "eta_scan") {
    return procData.eta_scan || null;
  }
  if (quantity.kind === "costheta_scan") {
    return procData.costheta_scan || null;
  }
  if (quantity.kind === "mass") {
    return procData.mass || null;
  }
  return null;
}

function currentQuantity() {
  const key = quantitySelect.value;
  return QUANTITIES.find(q => q.key === key);
}

function colorFor(method, process) {
  const k = method + "||" + process;
  if (selectedColors.has(k)) return selectedColors.get(k);
  const meta = DASHBOARD_DATA.process_meta[process] || {};
  const c = meta.color || AUTO_COLORS[autoColorIdx % AUTO_COLORS.length];
  autoColorIdx += 1;
  selectedColors.set(k, c);
  return c;
}

// --------------------------------------------------------------------------
// Presets: one-click combinations of methods/processes/colors mirroring the
// comparisons already made in joint_plots.py (jet multiplicity, clustering
// algo scan, detector/matching scan, energy-recovery on/off).
// --------------------------------------------------------------------------

function hslToHex(h, s, l) {
  s /= 100; l /= 100;
  const k = n => (n + h / 30) % 12;
  const a = s * Math.min(l, 1 - l);
  const f = n => l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
  const toHex = x => Math.round(255 * x).toString(16).padStart(2, "0");
  return `#${toHex(f(0))}${toHex(f(8))}${toHex(f(4))}`;
}

function antiKtRadius(methodName, recovery) {
  const re = recovery ? /^PF_E_recovery_AntiKtR(\d+)$/ : /^PF_AntiKtR(\d+)$/;
  const m = methodName.match(re);
  return m ? parseInt(m[1], 10) : null;
}

function representativeProcess() {
  return processNames.find(p => p.includes("qqqq")) || processNames[0];
}

function applyPreset(preset) {
  const methodSet = new Set(preset.methods);
  const processSet = new Set(preset.processes);
  methodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => {
    cb.checked = methodSet.has(cb.value);
  });
  processOptions.querySelectorAll("input[type=checkbox]").forEach(cb => {
    cb.checked = processSet.has(cb.value);
  });
  preset.methods.forEach(m => preset.processes.forEach(p => {
    const color = preset.colorOf(m, p);
    if (color) selectedColors.set(m + "||" + p, color);
  }));
  if (preset.quantity) quantitySelect.value = preset.quantity;
  redraw();
}

function buildPresets() {
  const presets = [];

  // 1. Jet multiplicity comparison: one method, every process, using the
  // existing process_config color families (teal=2 jets, magenta=4 jets, blue=6 jets).
  const durhamLike = methodNames.find(m => m === "PF_Durham") || methodNames[0];
  if (durhamLike) {
    presets.push({
      title: "2/4/6-jet processes",
      desc: `All processes, ${DASHBOARD_DATA.methods[durhamLike].label} — colored by jet multiplicity`,
      methods: [durhamLike],
      processes: processNames.slice(),
      colorOf: (m, p) => (DASHBOARD_DATA.process_meta[p] || {}).color,
      quantity: "energy:_all",
    });
  }

  // 2. Clustering algorithm scan: Durham + AntiKt radii (no energy recovery),
  // single representative process, AntiKt getting darker purple with radius.
  const akMethods = methodNames
    .filter(m => antiKtRadius(m, false) !== null)
    .sort((a, b) => antiKtRadius(a, false) - antiKtRadius(b, false));
  if (durhamLike && akMethods.length) {
    const proc = representativeProcess();
    const methods = [durhamLike, ...akMethods];
    const colorOf = (m, p) => {
      if (m === durhamLike) return "#1f77b4";
      const idx = akMethods.indexOf(m);
      const lightness = 70 - (idx / Math.max(1, akMethods.length - 1)) * 45;
      return hslToHex(270, 60, lightness);
    };
    presets.push({
      title: "Clustering algorithm scan",
      desc: `Durham vs. anti-kt radii, for ${(DASHBOARD_DATA.process_meta[proc] || {}).label || proc}`,
      methods, processes: [proc], colorOf, quantity: "energy:_all",
    });
  }

  // 3. Detector / matching comparison: PF vs Calo vs PF+IdealMatching.
  const detectorMethods = ["PF_Durham", "CaloJets_Durham", "PF_Durham_IdealMatching"].filter(
    m => methodNames.includes(m)
  );
  const detectorColors = {
    PF_Durham: "#1f77b4", CaloJets_Durham: "#2ca02c", PF_Durham_IdealMatching: "#ff7f0e",
  };
  if (detectorMethods.length > 1) {
    const proc = representativeProcess();
    presets.push({
      title: "Detector / matching comparison",
      desc: `PF jets vs. Calo jets vs. ideal matching, for ${(DASHBOARD_DATA.process_meta[proc] || {}).label || proc}`,
      methods: detectorMethods, processes: [proc],
      colorOf: (m) => detectorColors[m],
      quantity: "energy:_all",
    });
  }

  // 4. Energy recovery on/off: pair up AntiKt radii available both with and
  // without energy recovery; same hue per radius, recovery = saturated/dark,
  // no recovery = light, so pairs are visually grouped.
  const akNoRec = methodNames.filter(m => antiKtRadius(m, false) !== null);
  const akRec = methodNames.filter(m => antiKtRadius(m, true) !== null);
  const pairedRadii = akNoRec
    .map(m => antiKtRadius(m, false))
    .filter(r => akRec.some(m => antiKtRadius(m, true) === r))
    .sort((a, b) => a - b)
    .slice(0, 5);
  if (pairedRadii.length) {
    const proc = representativeProcess();
    const methods = [];
    const colorOf = (m) => {
      const rNoRec = antiKtRadius(m, false);
      const rRec = antiKtRadius(m, true);
      const r = rNoRec !== null ? rNoRec : rRec;
      const idx = pairedRadii.indexOf(r);
      const hue = (idx * 360) / Math.max(1, pairedRadii.length);
      return rNoRec !== null ? hslToHex(hue, 45, 75) : hslToHex(hue, 85, 40);
    };
    pairedRadii.forEach(r => {
      const noRec = methodNames.find(m => antiKtRadius(m, false) === r);
      const rec = methodNames.find(m => antiKtRadius(m, true) === r);
      if (noRec) methods.push(noRec);
      if (rec) methods.push(rec);
    });
    presets.push({
      title: "Energy recovery on/off",
      desc: `Anti-kt with vs. without energy recovery, by radius, for ${(DASHBOARD_DATA.process_meta[proc] || {}).label || proc}`,
      methods, processes: [proc], colorOf, quantity: "energy:_all",
    });
  }

  // 5. B-hadron content scan: single method, all processes, colored by the
  // existing LINE_STYLES convention (full = b-jets, dashed = mixed/gluons,
  // dotted = light-flavour), mirroring the presentation's "more B-hadron
  // content" row axis.
  if (durhamLike) {
    const flavourColor = { "-": "#8c2d04", "--": "#fd8d3c", ":": "#fdd0a2" };
    const flavourRank = { "-": 0, "--": 1, ":": 2 };
    presets.push({
      title: "B-hadron content scan",
      desc: `All processes, ${DASHBOARD_DATA.methods[durhamLike].label} — colored by B-hadron content (dark = b-jets, light = light-flavour)`,
      methods: [durhamLike],
      processes: processNames.slice().sort((a, b) => {
        const la = (DASHBOARD_DATA.process_meta[a] || {}).line_style || "-";
        const lb = (DASHBOARD_DATA.process_meta[b] || {}).line_style || "-";
        return (flavourRank[la] ?? 1) - (flavourRank[lb] ?? 1);
      }),
      colorOf: (m, p) => {
        const ls = (DASHBOARD_DATA.process_meta[p] || {}).line_style || "-";
        return flavourColor[ls] || "#999999";
      },
      quantity: "energy:_all",
    });
  }

  // 6. Clustering algorithm: Higgs mass comparison (Durham vs. anti-kt radii),
  // showing under/over-clustering directly in the mH peak (mirrors slides 19-21).
  if (durhamLike && akMethods.length) {
    const proc = representativeProcess();
    const methods = [durhamLike, ...akMethods];
    const colorOf = (m) => {
      if (m === durhamLike) return "#1f77b4";
      const idx = akMethods.indexOf(m);
      const lightness = 70 - (idx / Math.max(1, akMethods.length - 1)) * 45;
      return hslToHex(270, 60, lightness);
    };
    presets.push({
      title: "Higgs mass: clustering algorithm scan",
      desc: `Durham vs. anti-kt radii mH peak, for ${(DASHBOARD_DATA.process_meta[proc] || {}).label || proc}`,
      methods, processes: [proc], colorOf, quantity: "mass",
    });
  }

  // 7. Higgs mass: detector / matching comparison (mirrors slide 15).
  if (detectorMethods.length > 1) {
    const proc = representativeProcess();
    presets.push({
      title: "Higgs mass: detector comparison",
      desc: `PF jets vs. Calo jets vs. ideal matching mH peak, for ${(DASHBOARD_DATA.process_meta[proc] || {}).label || proc}`,
      methods: detectorMethods, processes: [proc],
      colorOf: (m) => detectorColors[m],
      quantity: "mass",
    });
  }

  const container = document.getElementById("presetList");
  container.innerHTML = "";
  presets.forEach(preset => {
    const btn = document.createElement("button");
    btn.className = "preset-btn";
    btn.innerHTML = `${preset.title}<span class="preset-desc">${preset.desc}</span>`;
    btn.addEventListener("click", () => applyPreset(preset));
    container.appendChild(btn);
  });
  if (presets.length === 0) {
    container.innerHTML = '<span style="font-size:11px;color:#999;">No presets available for this dataset.</span>';
  }
}

function rebuildSelectionList(combos) {
  const container = document.getElementById("selection-list");
  container.innerHTML = "";
  combos.forEach(({method, process}) => {
    const k = method + "||" + process;
    const row = document.createElement("div");
    row.className = "sel-row";
    const colorInput = document.createElement("input");
    colorInput.type = "color";
    colorInput.value = colorFor(method, process);
    colorInput.addEventListener("input", () => {
      selectedColors.set(k, colorInput.value);
      redraw();
    });
    const label = DASHBOARD_DATA.methods[method].label + " / " +
      ((DASHBOARD_DATA.process_meta[process] || {}).label || process);
    const span = document.createElement("span");
    span.textContent = label;
    span.title = label;
    row.appendChild(colorInput);
    row.appendChild(span);
    container.appendChild(row);
  });
}

let currentPointIndex = []; // parallel array: per main-plot trace point -> {method, process, quantity, binIdx}

function massHistTrace(rec, color, name, dash) {
  const edges = rec.edges || [];
  const y = rec.y || [];
  const centers = [];
  for (let i = 0; i < edges.length - 1; i++) centers.push(0.5 * (edges[i] + edges[i + 1]));
  return {
    x: centers, y: y, mode: "lines", type: "scatter",
    line: {color: color, dash: dash || "solid", shape: "hv"},
    name: name,
  };
}

function renderMassPlot(combos) {
  const quantity = currentQuantity();
  const traces = [];
  const showAllDefs = combos.length === 1;

  combos.forEach(({method, process}) => {
    const entry = getEntry(method, process, quantity);
    if (!entry) return;
    const baseColor = colorFor(method, process);
    const labelBase = DASHBOARD_DATA.methods[method].label + " / " +
      ((DASHBOARD_DATA.process_meta[process] || {}).label || process);
    const defsToShow = showAllDefs ? MASS_DEFINITIONS : MASS_DEFINITIONS.filter(d => d.key === "reco");

    defsToShow.forEach(def => {
      const rec = (entry.definitions || {})[def.key];
      if (!rec || !rec.edges || !rec.edges.length) return;
      const color = showAllDefs ? def.color : baseColor;
      const name = showAllDefs ? def.label : labelBase;
      traces.push(massHistTrace(rec, color, name));
    });

    if (entry.fit && entry.fit.fit_x && entry.fit.fit_y) {
      traces.push({
        x: entry.fit.fit_x, y: entry.fit.fit_y, mode: "lines", type: "scatter",
        name: (showAllDefs ? "reco" : labelBase) + " (Gaussian fit)",
        line: {color: showAllDefs ? "#000000" : baseColor, dash: "dash"},
        hoverinfo: "skip",
      });
    }
  });

  Plotly.react("mainPlot", traces, {
    xaxis: {title: "m_H [GeV]"},
    yaxis: {title: "Events"},
    margin: {t: 20},
    hovermode: "closest",
  }, {responsive: true});
}

function redraw() {
  const methods = getChecked(methodOptions);
  const processes = getChecked(processOptions);
  const quantity = currentQuantity();

  const combos = [];
  methods.forEach(m => processes.forEach(p => combos.push({method: m, process: p})));
  rebuildSelectionList(combos);

  currentPointIndex = [];

  const isMass = quantity.kind === "mass";
  document.getElementById("histHint").style.display = isMass ? "none" : "flex";
  document.getElementById("histPlot").style.display = isMass ? "none" : "block";
  document.getElementById("histSelectionList").style.display = isMass ? "none" : "block";
  if (isMass) {
    renderMassPlot(combos);
    return;
  }

  const traces = [];

  combos.forEach(({method, process}) => {
    const entry = getEntry(method, process, quantity);
    if (!entry) return;
    const color = colorFor(method, process);
    const labelBase = DASHBOARD_DATA.methods[method].label + " / " +
      ((DASHBOARD_DATA.process_meta[process] || {}).label || process);

    const yKey = entry.sigma_over_E !== undefined ? "sigma_over_E" : "sigma";
    const xs = entry.mid_points || [];
    const ys = entry[yKey] || [];

    const customdata = xs.map((_, i) => traces.length + "::" + i);

    traces.push({
      x: xs,
      y: ys,
      mode: "markers",
      type: "scatter",
      name: labelBase,
      marker: {color: color, size: 8},
      customdata: customdata,
    });
    currentPointIndex.push({method, process, quantity: quantity.key, entry});

    if (entry.fit_x && entry.fit_y) {
      traces.push({
        x: entry.fit_x,
        y: entry.fit_y,
        mode: "lines",
        type: "scatter",
        name: labelBase + " (fit)",
        line: {color: color, dash: "dash"},
        showlegend: false,
        hoverinfo: "skip",
      });
      currentPointIndex.push(null);
    }
  });

  const xAxisTitle = (quantity.kind === "eta_scan") ? "eta" :
    (quantity.kind === "costheta_scan") ? "cos(theta)" : "True jet energy [GeV]";
  const yAxisTitle = quantity.kind === "energy" || quantity.kind === "eta_scan" || quantity.kind === "costheta_scan"
    ? "sigma(E)/E" : "sigma";

  Plotly.react("mainPlot", traces, {
    xaxis: {title: xAxisTitle},
    yaxis: {title: yAxisTitle},
    margin: {t: 20},
    hovermode: "closest",
  }, {responsive: true});
}

// Selected histogram points, keyed by a stable id (survives redraws/re-coloring
// of the main plot, unlike trace/point indices). Map(key -> {info, binIdx, color}).
const selectedHistPoints = new Map();
const HIST_AUTO_COLORS = [
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
  "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
];

function histKey(info, binIdx) {
  return `${info.method}||${info.process}||${info.quantity}||${binIdx}`;
}

function attachClickHandler() {
  const gd = document.getElementById("mainPlot");
  gd.on("plotly_click", function(data) {
    const pt = data.points[0];
    const traceIdx = pt.curveNumber;
    const info = currentPointIndex[traceIdx];
    if (!info) return;
    const binIdx = pt.pointIndex;
    toggleHistPoint(info, binIdx);
  });
}

function toggleHistPoint(info, binIdx) {
  const key = histKey(info, binIdx);
  if (selectedHistPoints.has(key)) {
    selectedHistPoints.delete(key);
  } else {
    const color = HIST_AUTO_COLORS[selectedHistPoints.size % HIST_AUTO_COLORS.length];
    selectedHistPoints.set(key, {info, binIdx, color});
  }
  redrawHistPlot();
}

function rebuildHistSelectionList() {
  const container = document.getElementById("histSelectionList");
  container.innerHTML = "";
  selectedHistPoints.forEach((entry, key) => {
    const {info, binIdx, color} = entry;
    const rec = (info.entry.bins || [])[binIdx];
    const row = document.createElement("div");
    row.className = "hist-row";
    const swatch = document.createElement("div");
    swatch.className = "swatch";
    swatch.style.background = color;
    const label = `${DASHBOARD_DATA.methods[info.method].label} / ` +
      `${(DASHBOARD_DATA.process_meta[info.process] || {}).label || info.process} ` +
      `[${rec ? rec.lo : "?"}, ${rec ? rec.hi : "?"}]`;
    const span = document.createElement("span");
    span.textContent = label;
    span.title = label;
    const removeBtn = document.createElement("button");
    removeBtn.textContent = "✕";
    removeBtn.addEventListener("click", () => {
      selectedHistPoints.delete(key);
      redrawHistPlot();
    });
    row.appendChild(swatch);
    row.appendChild(span);
    row.appendChild(removeBtn);
    container.appendChild(row);
  });
}

function redrawHistPlot() {
  rebuildHistSelectionList();

  if (selectedHistPoints.size === 0) {
    Plotly.purge("histPlot");
    return;
  }

  const normalize = document.getElementById("normalizeHistCheckbox").checked;
  const traces = [];

  selectedHistPoints.forEach(({info, binIdx, color}) => {
    const rec = (info.entry.bins || [])[binIdx];
    if (!rec) return;
    const edges = rec.edges || [];
    let y = rec.y || [];
    const centers = [];
    const widths = [];
    for (let i = 0; i < edges.length - 1; i++) {
      centers.push(0.5 * (edges[i] + edges[i + 1]));
      widths.push(edges[i + 1] - edges[i]);
    }
    if (normalize) {
      const yMaxRaw = Math.max(1e-12, ...y);
      y = y.map(v => v / yMaxRaw);
    }

    const label = `${DASHBOARD_DATA.methods[info.method].label} / ` +
      `${(DASHBOARD_DATA.process_meta[info.process] || {}).label || info.process} ` +
      `[${rec.lo}, ${rec.hi}] GeV`;

    traces.push({
      x: centers,
      y: y,
      type: "bar",
      width: widths,
      marker: {color: color, opacity: 0.45},
      name: label,
      legendgroup: histKey(info, binIdx),
    });

    const yMax = normalize ? 1.0 : Math.max(1, ...y);
    const vline = (xval, name, dash) => ({
      x: [xval, xval], y: [0, yMax], mode: "lines", type: "scatter",
      name: `${name} (${label})`, line: {color: color, dash: dash || "solid"},
      legendgroup: histKey(info, binIdx), showlegend: false, hoverinfo: "skip",
    });

    if (rec.low !== undefined) traces.push(vline(rec.low, "low", "dot"));
    if (rec.high !== undefined) traces.push(vline(rec.high, "high", "dot"));
    if (rec.mpv !== undefined) traces.push(vline(rec.mpv, "MPV", "solid"));
  });

  Plotly.react("histPlot", traces, {
    barmode: "overlay",
    margin: {t: 20},
    showlegend: true,
    yaxis: {title: normalize ? "normalized" : "count"},
  }, {responsive: true});
}

// --------------------------------------------------------------------------
// Tabs
// --------------------------------------------------------------------------

document.querySelectorAll(".tab-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".tab-content").forEach(c => c.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById("tab-" + btn.dataset.tab).classList.add("active");
  });
});

// --------------------------------------------------------------------------
// Download raw JSON data (the page is self-contained, so this just re-offers
// the already-inlined DASHBOARD_DATA as a downloadable file).
// --------------------------------------------------------------------------

document.getElementById("downloadJsonBtn").addEventListener("click", () => {
  const blob = new Blob([JSON.stringify(DASHBOARD_DATA)], {type: "application/json"});
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "dashboard_data.json";
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
});

// --------------------------------------------------------------------------
// Statistics tab: fit coefficients, event counts, filter pass rates.
// --------------------------------------------------------------------------

function formatFitParams(model, popt) {
  if (model === "two_param" || model === "two_param_quadrature") {
    return `A=${popt[0].toFixed(3)} C=${popt[1].toFixed(3)}`;
  }
  if (model === "three_param" || model === "three_param_quadrature") {
    return `A=${popt[0].toFixed(3)} B=${popt[2].toFixed(3)} C=${popt[1].toFixed(3)}`;
  }
  if (model === "gaussian") {
    return `mu=${popt[0].toFixed(2)} sigma=${popt[1].toFixed(2)}`;
  }
  return popt.map(v => Number(v).toFixed(3)).join(", ");
}

function buildStatsTable(container, columns, rows) {
  if (!rows.length) {
    container.innerHTML = '<p style="font-size:12px;color:#999;">No data available.</p>';
    return;
  }
  let html = '<table class="stats-table"><thead><tr>';
  columns.forEach(c => { html += `<th>${c}</th>`; });
  html += "</tr></thead><tbody>";
  rows.forEach(row => {
    html += "<tr>";
    row.forEach(cell => {
      const isObj = cell !== null && typeof cell === "object";
      const cls = isObj && cell.best ? ' class="best"' : "";
      const text = isObj ? cell.text : cell;
      html += `<td${cls}>${text === undefined || text === null ? "" : text}</td>`;
    });
    html += "</tr>";
  });
  html += "</tbody></table>";
  container.innerHTML = html;
}

function renderFitCoefficientsTable() {
  const rows = [];
  methodNames.forEach(method => {
    const methodLabel = DASHBOARD_DATA.methods[method].label;
    processNames.forEach(process => {
      const procData = DASHBOARD_DATA.methods[method].processes[process];
      if (!procData) return;
      const procLabel = (DASHBOARD_DATA.process_meta[process] || {}).label || process;

      QUANTITIES.forEach(q => {
        if (q.kind === "mass") {
          const fit = procData.mass && procData.mass.fit;
          if (fit) {
            rows.push([methodLabel, procLabel, "Higgs mass peak", fit.model, formatFitParams(fit.model, fit.popt)]);
          }
          return;
        }
        const entry = getEntry(method, process, q);
        if (entry && entry.model && entry.popt) {
          rows.push([methodLabel, procLabel, q.label, entry.model, formatFitParams(entry.model, entry.popt)]);
        }
      });
    });
  });
  buildStatsTable(document.getElementById("fitCoeffTable"), ["Method", "Process", "Quantity", "Model", "Parameters"], rows);
}

function renderCountAndPassRateTables() {
  const stats = DASHBOARD_DATA.stats || {};
  const processes = Object.keys(stats).sort();
  const folderSet = new Set();
  processes.forEach(p => Object.keys(stats[p]).forEach(f => folderSet.add(f)));
  const folders = Array.from(folderSet).sort();

  const countRows = processes.map(p => {
    const row = [p];
    folders.forEach(f => {
      const m = stats[p][f];
      row.push(m ? `${Math.round(m.before)} / ${Math.round(m.after)}` : "");
    });
    return row;
  });
  buildStatsTable(document.getElementById("eventCountTable"), ["Process", ...folders], countRows);

  const rateRows = processes.map(p => {
    const values = folders.map(f => (stats[p][f] || {}).pass_rate);
    const numeric = values.filter(v => typeof v === "number");
    const best = numeric.length ? Math.max(...numeric) : null;
    const row = [p];
    folders.forEach((f, i) => {
      const v = values[i];
      if (typeof v !== "number") { row.push(""); return; }
      row.push({text: v.toFixed(3), best: best !== null && Math.abs(v - best) < 1e-12});
    });
    return row;
  });
  buildStatsTable(document.getElementById("passRateTable"), ["Process", ...folders], rateRows);
}

renderFitCoefficientsTable();
renderCountAndPassRateTables();

buildOptionList(methodOptions, methodNames, n => DASHBOARD_DATA.methods[n].label, redraw);
buildOptionList(processOptions, processNames, n => (DASHBOARD_DATA.process_meta[n] || {}).label || n, redraw);
buildPresets();
quantitySelect.addEventListener("change", redraw);

document.getElementById("clearHistBtn").addEventListener("click", () => {
  selectedHistPoints.clear();
  redrawHistPlot();
});
document.getElementById("normalizeHistCheckbox").addEventListener("change", redrawHistPlot);

// Sensible defaults: select first method/process so the dashboard isn't empty on load.
if (methodNames.length) methodOptions.querySelector("input[type=checkbox]").checked = true;
if (processNames.length) processOptions.querySelector("input[type=checkbox]").checked = true;

Plotly.newPlot("mainPlot", [], {margin: {t: 20}}).then(() => {
  attachClickHandler();
  redraw();
});
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data",
        type=str,
        required=True,
        help="Path to the dashboard_data.json produced by build_dashboard_data.py",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output HTML path (default: alongside --data, named dashboard.html)",
    )
    args = parser.parse_args()

    with open(args.data) as fh:
        data = json.load(fh)

    output_path = args.output or os.path.join(os.path.dirname(args.data), "dashboard.html")
    html = HTML_TEMPLATE.replace("__DASHBOARD_DATA_JSON__", json.dumps(data))

    with open(output_path, "w") as fh:
        fh.write(html)

    print("Saved interactive dashboard to:", output_path)


if __name__ == "__main__":
    main()
