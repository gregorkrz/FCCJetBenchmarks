"""Stage 4 of the resolution pipeline (no ROOT, no fitting — pure packaging).

Reads the consolidated `dashboard_data.json` produced by
`build_dashboard_data.py` and writes a single self-contained HTML file
(Plotly.js loaded from CDN, vanilla JS, no build step) with three tabs:

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
  can be compared side by side. Each drilled-down histogram is annotated with
  descriptive statistics (entry count N, mean, RMS, MPV, 68% band width), plus
  a one-line summary of how many histograms are plotted and their total entries.
- One-click presets (mirroring the comparisons in joint_plots.py and
  presentation.pdf) for common combinations: jet multiplicity (2/4/6 jets),
  clustering algorithm scan (Durham vs. anti-kt radii), detector/matching
  comparison (PF vs. Calo vs. ideal matching), energy recovery on/off,
  B-hadron content scan, and two Higgs-mass peak comparisons (clustering
  algorithm scan, detector comparison) — each with sensible default colors
  that can still be overridden afterwards.
- A "Download raw JSON data" button that re-offers the page's inlined
  DASHBOARD_DATA as a downloadable dashboard_data.json.

Grid tab: reproduces the 5x3 process matrix of joint_plots.py as a live grid -
one subplot per process, positioned by (B-hadron content x jet multiplicity),
with every selected clustering method / detector config overlaid inside each
cell. Grid presets mirror the reference presentation: detector/matching
comparison (PF vs. Calo vs. ideal matching), Higgs-mass clustering-algorithm
scan (Durham vs. anti-kt radii), and anti-kt energy recovery.

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
  .mini-btn { border: 1px solid #ccc; background: #f5f5f5; border-radius: 3px; padding: 2px 8px;
              font-size: 11px; cursor: pointer; margin: 0 4px 4px 0; }
  .mini-btn:hover { background: #eee; }
  /* Grid tab */
  #grid-layout { display: flex; gap: 16px; align-items: flex-start; }
  #grid-controls { flex: 0 0 300px; background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 12px; }
  #grid-plots-wrap { flex: 1; min-width: 0; }
  #grid-topaxis { text-align: center; font-size: 12px; font-weight: bold; color: #555; padding: 2px 0 6px 0; }
  #grid-body { display: flex; align-items: stretch; }
  #grid-leftaxis { flex: 0 0 22px; display: flex; align-items: center; justify-content: center; }
  #grid-leftaxis span { writing-mode: vertical-rl; transform: rotate(180deg);
                        font-size: 12px; font-weight: bold; color: #555; white-space: nowrap; }
  #gridPlots { flex: 1; min-width: 0; display: grid; gap: 6px; }
  .grid-cell { border: 1px solid #eee; border-radius: 4px; background: #fff; overflow: hidden; }
  .grid-cell .grid-cell-title { font-size: 11px; font-weight: bold; text-align: center;
                                padding: 3px 2px 0 2px; color: #333; white-space: nowrap;
                                overflow: hidden; text-overflow: ellipsis; }
  .grid-cell .grid-cell-plot { width: 100%; height: 190px; }
  .grid-cell.empty { border: 1px dashed #eee; background: #fafafa; }
  #grid-legend { display: flex; flex-wrap: wrap; gap: 10px; font-size: 12px; padding: 8px 0 0 0; }
  #grid-legend .leg-item { display: flex; align-items: center; gap: 5px; }
  #grid-legend .leg-swatch { width: 14px; height: 3px; border-radius: 2px; }
  .hist-stats { font-size: 10px; color: #666; margin-left: 18px; }
  #histStatsSummary { font-size: 11px; color: #555; padding: 4px 0; }
  /* Explanatory text for first-time users */
  #intro { background: #eef4fb; border: 1px solid #cfe0f3; border-radius: 6px;
           padding: 10px 14px; margin: 0 0 12px 0; font-size: 13px; color: #234; line-height: 1.5; }
  #intro h2 { margin: 0 0 4px 0; font-size: 15px; }
  #intro ul { margin: 6px 0 0 0; padding-left: 20px; }
  #intro li { margin: 2px 0; }
  #intro details { margin-top: 6px; }
  #intro summary { cursor: pointer; font-weight: bold; color: #1a4a80; }
  .tab-help { background: #f7f9fc; border-left: 3px solid #9cc0ea; border-radius: 3px;
              padding: 7px 11px; margin: 0 0 12px 0; font-size: 12.5px; color: #345; line-height: 1.45; }
  .res-warning { background: #fff6e6; border: 1px solid #f0d9a8; border-radius: 5px;
                 padding: 7px 11px; margin: 6px 0; font-size: 12px; color: #6b4e12; line-height: 1.45; }
  .res-warning b { color: #7a4e00; }
  #coeffModelDesc { background: #f4f7f4; border: 1px solid #d6e2d6; border-radius: 6px;
                    padding: 10px 14px; margin: 0 0 12px 0; font-size: 12.5px; color: #243; line-height: 1.5; }
  #coeffModelDesc h3 { margin: 0 0 5px 0; font-size: 14px; }
  #coeffModelDesc p { margin: 6px 0; }
  #coeffModelDesc .formula { font-family: ui-monospace, Menlo, Consolas, monospace; background: #e8efe8;
                             padding: 1px 5px; border-radius: 3px; }
  #coeffModelDesc .metrics { font-size: 11.5px; color: #567; margin-top: 6px; }
</style>
</head>
<body>
<div id="header-row">
  <h1>FCC Jet Resolution Dashboard</h1>
  <button id="downloadJsonBtn">Download raw JSON data</button>
</div>
<div id="intro">
  <h2>What is this?</h2>
  Interactive browser for <b>jet reconstruction performance</b> in FCC-ee ZH events at
  &radic;s = 240 GeV: jet <b>energy</b> and <b>angular resolution</b> vs. true jet energy, and
  reconstructed <b>Higgs-mass</b> peaks &mdash; compared across jet clustering methods
  (Particle-Flow / Calo, Durham / anti-k<sub>T</sub> at several radii, with/without energy recovery)
  and physics processes. Everything runs locally in your browser; nothing is uploaded.
  <details>
    <summary>How to use it (click to expand)</summary>
    <ul>
      <li><b>Three tabs</b> (top): <b>Explorer</b> = one overlay plot you build yourself;
          <b>Grid</b> = all processes laid out in the 2/4/6-jet &times; b-content matrix;
          <b>Statistics</b> = fit-coefficient and event-count tables.</li>
      <li><b>Presets</b> (left) are one-click starting points. Or pick a <b>Quantity</b>, then
          tick <b>Methods</b> and <b>Processes</b> to overlay.</li>
      <li><b>Click any data point</b> on a resolution plot to show the underlying histogram
          for that energy bin below; click more points to compare. Click a swatch to recolour.</li>
      <li><b>Fit model</b> dropdown switches the fitted curve (e.g. the classic
          A/&radic;E+B+C/E vs. the better log(E) forms found by our fit study). The
          <b>Statistics</b> tab lists every fit's coefficients and &chi;&sup2;/ndf.</li>
      <li>Your current view is saved in the page URL &mdash; copy it to share the exact plot.</li>
    </ul>
  </details>
</div>

<div id="tabBar">
  <button class="tab-btn active" data-tab="explorer">Explorer</button>
  <button class="tab-btn" data-tab="grid">Grid</button>
  <button class="tab-btn" data-tab="coeff">Coefficients</button>
  <button class="tab-btn" data-tab="statistics">Statistics</button>
</div>

<div id="tab-explorer" class="tab-content active">
<div class="tab-help">
  <b>Explorer.</b> Choose a quantity and tick methods/processes on the left to overlay them.
  <b>Click a marker</b> on the plot to reveal that energy bin's histogram below (click several to
  compare). Toggle <b>Show fit curve(s)</b>, <b>Show error bars</b>, or switch the <b>Fit model</b>.
</div>
<div id="layout">
  <div id="controls">
    <fieldset>
      <legend>Presets</legend>
      <div id="presetList"></div>
    </fieldset>
    <fieldset>
      <legend>Quantity</legend>
      <select id="quantitySelect" style="width: 100%;"></select>
      <label class="opt" style="margin-top: 6px;">
        <input type="checkbox" id="showFitCheckbox" checked> Show fit curve(s)
      </label>
      <label class="opt">
        <input type="checkbox" id="showErrorBarsCheckbox"> Show error bars (sigma/E &divide; sqrt(2N))
      </label>
      <label class="opt">
        Fit model:
        <select id="fitModelSelect" style="width: 100%;"></select>
      </label>
      <label class="opt">
        <input type="checkbox" id="showFullHistCheckbox"> Show full-resolution histograms (loads on demand)
      </label>
      <span id="fullHistStatus" style="font-size:11px;color:#a00;"></span>
    </fieldset>
    <fieldset>
      <legend>Methods</legend>
      <div>
        <button type="button" class="mini-btn" id="methodsAllBtn">All</button>
        <button type="button" class="mini-btn" id="methodsNoneBtn">None</button>
      </div>
      <div class="scroll-box" id="methodOptions"></div>
    </fieldset>
    <fieldset>
      <legend>Processes</legend>
      <div>
        <button type="button" class="mini-btn" id="processesAllBtn">All</button>
        <button type="button" class="mini-btn" id="processesNoneBtn">None</button>
        <button type="button" class="mini-btn" id="processes2jBtn">2-jet</button>
        <button type="button" class="mini-btn" id="processes4jBtn">4-jet</button>
        <button type="button" class="mini-btn" id="processes6jBtn">6-jet</button>
      </div>
      <div class="scroll-box" id="processOptions"></div>
    </fieldset>
    <fieldset>
      <legend>Selected curves (click swatch to recolor)</legend>
      <button type="button" class="mini-btn" id="autoColorsBtn">Auto colors</button>
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
    <div class="res-warning">
      &#9888; The histograms shown here are <b>downsampled previews</b> (rebinned for a small,
      fast page) &mdash; fine for eyeballing shapes, but not the full binning. For the exact
      histogram, tick <b>&ldquo;Show full-resolution histograms&rdquo;</b> in the Quantity box on
      the left; it fetches the full-binning data on demand for the points you've selected.
      All resolution/fit numbers are always computed from the full-resolution data, never the preview.
    </div>
    <div id="histStatsSummary"></div>
    <div id="histSelectionList"></div>
  </div>
</div>
</div>

<div id="tab-grid" class="tab-content">
<div class="tab-help">
  <b>Grid.</b> Every process shown at once, laid out in the physics matrix: columns are jet
  multiplicity (2 &rarr; 4 &rarr; 6 final-state jets), rows increase in b-hadron content. Pick a
  <b>quantity</b> and a <b>preset</b> (or tick methods) to compare clustering algorithms across all
  processes at a glance.
</div>
<div id="grid-layout">
  <div id="grid-controls">
    <fieldset>
      <legend>Grid presets</legend>
      <div id="gridPresetList"></div>
    </fieldset>
    <fieldset>
      <legend>Quantity</legend>
      <select id="gridQuantitySelect" style="width: 100%;"></select>
      <label class="opt" style="margin-top: 6px;">
        <input type="checkbox" id="gridShowFitCheckbox" checked> Show fit curve(s)
      </label>
      <label class="opt">
        Fit model:
        <select id="gridFitModelSelect" style="width: 100%;"></select>
      </label>
    </fieldset>
    <fieldset>
      <legend>Methods (overlaid in each cell)</legend>
      <div>
        <button type="button" class="mini-btn" id="gridMethodsAllBtn">All</button>
        <button type="button" class="mini-btn" id="gridMethodsNoneBtn">None</button>
        <button type="button" class="mini-btn" id="gridMethodsAutoColorBtn">Auto colors</button>
      </div>
      <div class="scroll-box" id="gridMethodOptions"></div>
    </fieldset>
    <fieldset>
      <legend>Method colors (click swatch to recolor)</legend>
      <div class="scroll-box" id="gridColorList"></div>
    </fieldset>
  </div>
  <div id="grid-plots-wrap">
    <div id="grid-topaxis">Higher number of final-state jets →</div>
    <div id="grid-body">
      <div id="grid-leftaxis"><span>More B-hadron content →</span></div>
      <div id="gridPlots"></div>
    </div>
    <div id="grid-legend"></div>
  </div>
</div>
</div>

<div id="tab-coeff" class="tab-content">
<div class="tab-help">
  <b>Coefficients.</b> How each <b>fitted parameter</b> of the jet-energy-resolution curve trends
  across configurations. Pick a <b>fit model</b> and an <b>x-axis</b> (anti-k<sub>T</sub> radius,
  number of final-state jets, or clustering method); you get one panel per coefficient
  (A, B, C, D&hellip;), with error bars from the fit covariance. Tick which methods/processes to
  include on the left.
</div>
<div id="coeff-layout" style="display:flex; gap:16px; align-items:flex-start;">
  <div id="coeff-controls" style="flex:0 0 300px; background:#fff; border:1px solid #ddd; border-radius:6px; padding:12px;">
    <fieldset>
      <legend>Fit model</legend>
      <select id="coeffFitModelSelect" style="width:100%;"></select>
    </fieldset>
    <fieldset>
      <legend>X-axis</legend>
      <select id="coeffXaxisSelect" style="width:100%;">
        <option value="radius">anti-kT radius R</option>
        <option value="njets">number of final-state jets</option>
        <option value="method">clustering method</option>
      </select>
      <div id="coeffXaxisHint" style="font-size:11px; color:#777; margin-top:5px;"></div>
    </fieldset>
    <fieldset>
      <legend>Methods</legend>
      <div>
        <button type="button" class="mini-btn" id="coeffMethodsAllBtn">All</button>
        <button type="button" class="mini-btn" id="coeffMethodsNoneBtn">None</button>
      </div>
      <div class="scroll-box" id="coeffMethodOptions"></div>
    </fieldset>
    <fieldset>
      <legend>Processes</legend>
      <div>
        <button type="button" class="mini-btn" id="coeffProcessesAllBtn">All</button>
        <button type="button" class="mini-btn" id="coeffProcessesNoneBtn">None</button>
        <button type="button" class="mini-btn" id="coeffProcesses2jBtn">2-jet</button>
        <button type="button" class="mini-btn" id="coeffProcesses4jBtn">4-jet</button>
        <button type="button" class="mini-btn" id="coeffProcesses6jBtn">6-jet</button>
      </div>
      <div class="scroll-box" id="coeffProcessOptions"></div>
    </fieldset>
  </div>
  <div id="coeff-plots-wrap" style="flex:1; min-width:0;">
    <div id="coeffModelDesc"></div>
    <div id="coeffPlots" style="display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:12px;"></div>
    <div id="coeffLegend" style="display:flex; flex-wrap:wrap; gap:12px; font-size:12px; padding:10px 0 0 0;"></div>
    <div id="coeffMsg" style="font-size:12px; color:#a00; padding:6px 0;"></div>
  </div>
</div>
</div>

<div id="tab-statistics" class="tab-content">
  <div class="tab-help">
    <b>Statistics.</b> Every fitted curve's coefficients (with uncertainties) and &chi;&sup2;/ndf,
    for all methods, processes and fit models, plus raw event counts and filter pass rates.
    Note: because the per-bin statistics are enormous (millions of jets), the statistical error
    bars are tiny, so &chi;&sup2;/ndf values are large across the board &mdash; use them to
    <em>rank</em> fits relative to each other, not as an absolute goodness-of-fit.
  </div>
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

// Energy-JER fit models offered in the fit-model dropdowns. `key` must match a
// key in each entry's `fits` object (produced by build_dashboard_data.py).
// Ordered best-first by leave-one-out CV (see fit_trials/README.md); the classic
// forms are kept at the end for reference. Both Explorer and Grid dropdowns are
// populated from this one list so they never drift apart.
const FIT_MODELS = [
  {key: "logE_4param_weighted", label: "log(E) 4-param, weighted (CV best)"},
  {key: "logE_4param",          label: "log(E) 4-param: A/√E+B+C/E+D·logE"},
  {key: "invE15_noise",         label: "4-param + e/E^1.5"},
  {key: "invE2_noise",          label: "4-param + e/E^2"},
  {key: "invE2_noise_weighted", label: "4-param + e/E^2, weighted"},
  {key: "logE_invE2_noStoch",   label: "no-√E: B+C/E+e/E^2+D·logE"},
  {key: "logE_and_invE2",       label: "5-param log+e/E^2 (best in-sample)"},
  {key: "logE_and_invE2_weighted", label: "5-param log+e/E^2, weighted"},
  {key: "three_param",          label: "3-param (classic): A/√E+B+C/E"},
  {key: "three_param_weighted", label: "3-param, weighted"},
  {key: "two_param",            label: "2-param: A/√E+B"},
  {key: "two_param_weighted",   label: "2-param, weighted"},
  {key: "staged_const_first",   label: "3-param staged: constant pinned to high-E tail"},
];

// Dashboard fit keys whose stored popt follows another model's parameter
// convention - staged fits reuse the classic three_param algebra and only differ
// in how the parameters are obtained, so they share its labels and formatting.
const FIT_BASE_ALIAS = {
  staged_const_first: "three_param",
};
function fitBaseModelOf(model) {
  const base = model.replace(/_weighted$/, "");
  return FIT_BASE_ALIAS[base] || base;
}

function populateFitModelSelect(sel, defaultKey) {
  FIT_MODELS.forEach(fm => {
    const opt = document.createElement("option");
    opt.value = fm.key;
    opt.textContent = fm.label;
    sel.appendChild(opt);
  });
  sel.value = defaultKey;
}
populateFitModelSelect(document.getElementById("fitModelSelect"), "logE_4param_weighted");
populateFitModelSelect(document.getElementById("gridFitModelSelect"), "logE_4param_weighted");

const methodOptions = document.getElementById("methodOptions");
const processOptions = document.getElementById("processOptions");

function getEntry(method, process, quantity) {
  const methodData = DASHBOARD_DATA.methods && DASHBOARD_DATA.methods[method];
  const procData = methodData && methodData.processes && methodData.processes[process];
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

// --------------------------------------------------------------------------
// Full-resolution histograms: each light bin/definition record carries its
// own "full_url" (written by build_dashboard_data.py into a `full_hist/`
// folder next to dashboard_data.json / dashboard.html), one small JSON file
// per histogram. Only the exact histogram currently being rendered is fetched
// - never the whole full-resolution dataset - and only once the user asks to
// see non-downsampled histograms.
// --------------------------------------------------------------------------

const fullRecordCache = new Map(); // full_url -> {edges, y}
const fullRecordPromises = new Map(); // full_url -> in-flight promise
const pendingFullUrls = new Set();

function updateFullHistStatus() {
  const status = document.getElementById("fullHistStatus");
  status.textContent = pendingFullUrls.size
    ? `Loading ${pendingFullUrls.size} full-resolution histogram(s)...`
    : "";
}

function showFullHistEnabled() {
  return document.getElementById("showFullHistCheckbox").checked;
}

// Given a light bin/definition record (with "edges"/"y" and, if a
// full-resolution version exists, "full_url"), return the full-resolution
// version if the user asked for it and it's already cached; otherwise kick
// off (at most one) fetch for that exact histogram and return the light
// version in the meantime, redrawing once the fetch resolves.
function resolveHistRecord(lightRec) {
  if (!showFullHistEnabled() || !lightRec || !lightRec.full_url) return lightRec;
  const url = lightRec.full_url;
  if (fullRecordCache.has(url)) return fullRecordCache.get(url);
  if (!fullRecordPromises.has(url)) {
    pendingFullUrls.add(url);
    updateFullHistStatus();
    const promise = fetch(url)
      .then(r => { if (!r.ok) throw new Error("HTTP " + r.status); return r.json(); })
      .then(data => {
        fullRecordCache.set(url, data);
        fullRecordPromises.delete(url);
        pendingFullUrls.delete(url);
        updateFullHistStatus();
        redraw();
        redrawHistPlot();
      })
      .catch(err => {
        fullRecordPromises.delete(url);
        pendingFullUrls.delete(url);
        updateFullHistStatus();
        document.getElementById("fullHistStatus").textContent =
          `Could not load ${url} (${err.message}). Make sure full_hist/ is next to dashboard.html ` +
          "and both are served over http(s), not file://.";
      });
    fullRecordPromises.set(url, promise);
  }
  return lightRec;
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
    const flavourColor = { "-": "#08519c", "--": "#4292c6", ":": "#c6dbef" };
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

function showFitEnabled() {
  return document.getElementById("showFitCheckbox").checked;
}

function showErrorBarsEnabled() {
  return document.getElementById("showErrorBarsCheckbox").checked;
}

function selectedFitModel() {
  return document.getElementById("fitModelSelect").value;
}

// Resolve which fit to draw for an energy/angle entry, given the selected
// fit-model dropdown. Prefers entry.fits[model] (both models pre-fit in
// resolution_plots.py); falls back to the entry's own single model/popt for
// older-shaped data that only has one fit.
function resolveFit(entry) {
  const model = selectedFitModel();
  if (entry.fits && entry.fits[model]) return entry.fits[model];
  if (entry.model === model && entry.fit_x && entry.fit_y) {
    return {popt: entry.popt, fit_x: entry.fit_x, fit_y: entry.fit_y};
  }
  return null;
}

function renderMassPlot(combos) {
  const quantity = currentQuantity();
  const traces = [];
  const showAllDefs = combos.length === 1;
  const showFit = showFitEnabled();

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
      const histRec = resolveHistRecord(rec);
      const color = showAllDefs ? def.color : baseColor;
      const name = showAllDefs ? def.label : labelBase;
      traces.push(massHistTrace(histRec, color, name));
    });

    if (showFit && entry.fit && entry.fit.fit_x && entry.fit.fit_y) {
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
    saveStateToHash();
    return;
  }

  const traces = [];
  const showFit = showFitEnabled();

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

    const trace = {
      x: xs,
      y: ys,
      mode: "markers",
      type: "scatter",
      name: labelBase,
      marker: {color: color, size: 8},
      customdata: customdata,
    };
    // Error bars (energy JER only): sigma_over_E_err is baked in by
    // build_dashboard_data.py as delta(sigma/E) = (sigma/E)/sqrt(2N). Gated by
    // the checkbox; other quantities simply have no error array to show.
    if (showErrorBarsEnabled() && Array.isArray(entry.sigma_over_E_err)) {
      trace.error_y = {
        type: "data",
        array: entry.sigma_over_E_err,
        visible: true,
        thickness: 1,
        width: 3,
        color: color,
      };
    }
    traces.push(trace);
    currentPointIndex.push({method, process, quantity: quantity.key, entry});

    const fit = showFit ? resolveFit(entry) : null;
    if (fit) {
      traces.push({
        x: fit.fit_x,
        y: fit.fit_y,
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

  saveStateToHash();
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

// Descriptive statistics for a single drilled-down histogram: the number of
// jets/entries that went into the bin (rec.n_jets, carried straight through
// from resolution_plots.py), plus the histogram's own mean / RMS / integral
// computed from its (bin-center, weight) contents, and the pre-computed 68%
// band (low/high) and most-probable value (MPV). Uses the light record for the
// metadata; falls back to the full-resolution shape when it's loaded.
function histStats(rec) {
  const histRec = resolveHistRecord(rec);
  const edges = histRec.edges || [];
  const y = histRec.y || [];
  let sumW = 0, sumWX = 0, sumWXX = 0, peakY = -Infinity, peakX = null;
  for (let i = 0; i < y.length && i < edges.length - 1; i++) {
    const c = 0.5 * (edges[i] + edges[i + 1]);
    const w = y[i];
    sumW += w; sumWX += w * c; sumWXX += w * c * c;
    if (w > peakY) { peakY = w; peakX = c; }
  }
  const mean = sumW > 0 ? sumWX / sumW : null;
  const variance = sumW > 0 ? Math.max(0, sumWXX / sumW - mean * mean) : null;
  const rms = variance !== null ? Math.sqrt(variance) : null;
  const width = (rec.low !== undefined && rec.high !== undefined) ? rec.high - rec.low : null;
  return {
    nEntries: rec.n_jets,
    mean, rms, peakX,
    integral: sumW,
    mpv: rec.mpv,
    low: rec.low,
    high: rec.high,
    bandWidth: width,
  };
}

function fmtNum(v, digits) {
  if (v === undefined || v === null || Number.isNaN(v)) return "-";
  return Number(v).toFixed(digits === undefined ? 3 : digits);
}

function histStatsText(rec) {
  const s = histStats(rec);
  const parts = [];
  if (s.nEntries !== undefined && s.nEntries !== null) {
    parts.push(`N=${Math.round(s.nEntries).toLocaleString()}`);
  }
  if (s.mean !== null) parts.push(`mean=${fmtNum(s.mean)}`);
  if (s.rms !== null) parts.push(`RMS=${fmtNum(s.rms)}`);
  if (s.mpv !== undefined) parts.push(`MPV=${fmtNum(s.mpv)}`);
  if (s.bandWidth !== null) parts.push(`68%-width=${fmtNum(s.bandWidth)}`);
  return parts.join(" · ");
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
    if (rec) {
      const stats = document.createElement("span");
      stats.className = "hist-stats";
      stats.textContent = histStatsText(rec);
      span.appendChild(document.createElement("br"));
      span.appendChild(stats);
    }
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

// One-line summary across all currently-drilled-down histograms: how many are
// shown and the total number of jets/entries behind them.
function updateHistStatsSummary() {
  const el = document.getElementById("histStatsSummary");
  if (!el) return;
  if (selectedHistPoints.size === 0) { el.textContent = ""; return; }
  let totalEntries = 0, counted = 0;
  selectedHistPoints.forEach(({info, binIdx}) => {
    const rec = (info.entry.bins || [])[binIdx];
    if (rec && rec.n_jets !== undefined && rec.n_jets !== null) {
      totalEntries += rec.n_jets; counted += 1;
    }
  });
  const n = selectedHistPoints.size;
  let text = `${n} histogram${n === 1 ? "" : "s"} plotted`;
  if (counted) text += ` · ${Math.round(totalEntries).toLocaleString()} jets/entries total`;
  el.textContent = text;
}

function redrawHistPlot() {
  rebuildHistSelectionList();
  updateHistStatsSummary();

  if (selectedHistPoints.size === 0) {
    Plotly.purge("histPlot");
    saveStateToHash();
    return;
  }

  const normalize = document.getElementById("normalizeHistCheckbox").checked;
  const traces = [];

  selectedHistPoints.forEach(({info, binIdx, color}) => {
    const rec = (info.entry.bins || [])[binIdx];
    if (!rec) return;
    // lo/hi/low/high/mpv metadata only exists on the light record - the full
    // version (when loaded) only carries edges/y, so it's resolved separately
    // and only used for the histogram shape itself.
    const histRec = resolveHistRecord(rec);
    const edges = histRec.edges || [];
    let y = histRec.y || [];
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

  saveStateToHash();
}

// --------------------------------------------------------------------------
// URL state persistence: the full UI state (selected methods/processes,
// quantity, fit controls, colors, drilled-down histogram points, active tab)
// is base64url-encoded into the URL hash after every change, and restored
// from it on load - so a specific dashboard view can be bookmarked or shared
// as a link. Uses history.replaceState (not location.hash=...) so writing
// the state doesn't spam browser history or fire a "hashchange" event for
// our own writes.
// --------------------------------------------------------------------------

function encodeState(obj) {
  const json = JSON.stringify(obj);
  const b64 = btoa(unescape(encodeURIComponent(json)));
  // base64url, no padding: a "clean" alphabet (A-Za-z0-9-_) that never needs
  // percent-encoding inside a URL hash.
  return b64.replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
}

function decodeState(str) {
  if (!str) return null;
  try {
    let b64 = str.replace(/-/g, "+").replace(/_/g, "/");
    while (b64.length % 4) b64 += "=";
    const json = decodeURIComponent(escape(atob(b64)));
    return JSON.parse(json);
  } catch (e) {
    console.warn("Could not decode dashboard state from URL:", e);
    return null;
  }
}

function captureState() {
  const activeTabBtn = document.querySelector(".tab-btn.active");
  return {
    tab: activeTabBtn ? activeTabBtn.dataset.tab : "explorer",
    quantity: quantitySelect.value,
    fitModel: document.getElementById("fitModelSelect").value,
    showFit: document.getElementById("showFitCheckbox").checked,
    showErrorBars: document.getElementById("showErrorBarsCheckbox").checked,
    showFullHist: document.getElementById("showFullHistCheckbox").checked,
    normalizeHist: document.getElementById("normalizeHistCheckbox").checked,
    methods: getChecked(methodOptions),
    processes: getChecked(processOptions),
    colors: Object.fromEntries(selectedColors),
    histPoints: Array.from(selectedHistPoints.values()).map(({info, binIdx}) => ({
      method: info.method, process: info.process, quantity: info.quantity, binIdx,
    })),
    grid: {
      quantity: gridQuantitySelect.value,
      showFit: document.getElementById("gridShowFitCheckbox").checked,
      fitModel: document.getElementById("gridFitModelSelect").value,
      methods: getChecked(gridMethodOptions),
      colors: Object.fromEntries(gridMethodColors),
    },
    coeff: {
      fitModel: coeffFitModelSelect.value,
      xaxis: coeffXaxisSelect.value,
      methods: getChecked(coeffMethodOptions),
      processes: getChecked(coeffProcessOptions),
    },
  };
}

function applyState(state) {
  if (!state) return;
  if (state.quantity) quantitySelect.value = state.quantity;
  if (state.fitModel) document.getElementById("fitModelSelect").value = state.fitModel;
  if (typeof state.showFit === "boolean") document.getElementById("showFitCheckbox").checked = state.showFit;
  if (typeof state.showErrorBars === "boolean") document.getElementById("showErrorBarsCheckbox").checked = state.showErrorBars;
  if (typeof state.normalizeHist === "boolean") {
    document.getElementById("normalizeHistCheckbox").checked = state.normalizeHist;
  }

  const methodSet = new Set(state.methods || []);
  const processSet = new Set(state.processes || []);
  methodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = methodSet.has(cb.value); });
  processOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = processSet.has(cb.value); });

  selectedColors.clear();
  Object.entries(state.colors || {}).forEach(([k, v]) => selectedColors.set(k, v));

  selectedHistPoints.clear();
  (state.histPoints || []).forEach(({method, process, quantity, binIdx}) => {
    const q = QUANTITIES.find(qq => qq.key === quantity);
    if (!q) return;
    const entry = getEntry(method, process, q);
    if (!entry) return;
    const info = {method, process, quantity, entry};
    const key = histKey(info, binIdx);
    const color = HIST_AUTO_COLORS[selectedHistPoints.size % HIST_AUTO_COLORS.length];
    selectedHistPoints.set(key, {info, binIdx, color});
  });

  const g = state.grid || {};
  if (g.quantity) gridQuantitySelect.value = g.quantity;
  if (g.fitModel) document.getElementById("gridFitModelSelect").value = g.fitModel;
  if (typeof g.showFit === "boolean") document.getElementById("gridShowFitCheckbox").checked = g.showFit;
  const gridMethodSet = new Set(g.methods || []);
  gridMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => {
    cb.checked = gridMethodSet.has(cb.value);
  });
  gridMethodColors.clear();
  Object.entries(g.colors || {}).forEach(([k, v]) => gridMethodColors.set(k, v));

  const co = state.coeff || {};
  if (co.fitModel) coeffFitModelSelect.value = co.fitModel;
  if (co.xaxis) coeffXaxisSelect.value = co.xaxis;
  if (co.methods) {
    const s = new Set(co.methods);
    coeffMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = s.has(cb.value); });
  }
  if (co.processes) {
    const s = new Set(co.processes);
    coeffProcessOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = s.has(cb.value); });
  }

  if (state.tab === "statistics" || state.tab === "grid" || state.tab === "coeff") {
    document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".tab-content").forEach(c => c.classList.remove("active"));
    const btn = Array.from(document.querySelectorAll(".tab-btn")).find(b => b.dataset.tab === state.tab);
    if (btn) btn.classList.add("active");
    document.getElementById("tab-" + state.tab).classList.add("active");
  }

  if (typeof state.showFullHist === "boolean") {
    document.getElementById("showFullHistCheckbox").checked = state.showFullHist;
  }
  redraw();
  redrawHistPlot();
  if (state.tab === "grid") redrawGrid();
  if (state.tab === "coeff") redrawCoeff();
}

function saveStateToHash() {
  const encoded = encodeState(captureState());
  history.replaceState(null, "", "#" + encoded);
}

// --------------------------------------------------------------------------
// Grid tab: reproduces the 5x3 process matrix of joint_plots.py as a live,
// interactive grid. One subplot per process, positioned by its
// (grid_row = B-hadron content, grid_col = jet multiplicity) coordinates, with
// every selected method overlaid inside each cell - so the whole method
// comparison (e.g. PF vs. Calo vs. ideal matching, or a clustering-algorithm /
// energy-recovery scan of the Higgs mass) can be read across all processes at
// once. Presets mirror slides 12/16/19/21 of the reference presentation.
// --------------------------------------------------------------------------

const gridQuantitySelect = document.getElementById("gridQuantitySelect");
QUANTITIES.forEach(q => {
  const opt = document.createElement("option");
  opt.value = q.key;
  opt.textContent = q.label;
  gridQuantitySelect.appendChild(opt);
});

const gridMethodOptions = document.getElementById("gridMethodOptions");

// Method-keyed colors (the grid overlays methods, unlike the Explorer which
// keys colors by method||process). Map(method -> hex color).
const gridMethodColors = new Map();
let gridAutoColorIdx = 0;

function gridColorForMethod(method) {
  if (gridMethodColors.has(method)) return gridMethodColors.get(method);
  const c = AUTO_COLORS[gridAutoColorIdx % AUTO_COLORS.length];
  gridAutoColorIdx += 1;
  gridMethodColors.set(method, c);
  return c;
}

function gridCurrentQuantity() {
  return QUANTITIES.find(q => q.key === gridQuantitySelect.value);
}

function gridSelectedMethods() {
  return getChecked(gridMethodOptions);
}

// Processes that have a grid position, keyed "row,col" -> process name, plus
// the matrix extent. Fixed for a given dataset (driven by process_config.py's
// PROCESS_TO_ROW_COL, exported into process_meta by build_dashboard_data.py).
function gridLayout() {
  const cellOf = {};
  let nRows = 0, nCols = 0;
  processNames.forEach(p => {
    const meta = DASHBOARD_DATA.process_meta[p] || {};
    if (meta.grid_row === null || meta.grid_row === undefined) return;
    if (meta.grid_col === null || meta.grid_col === undefined) return;
    cellOf[meta.grid_row + "," + meta.grid_col] = p;
    nRows = Math.max(nRows, meta.grid_row + 1);
    nCols = Math.max(nCols, meta.grid_col + 1);
  });
  return {cellOf, nRows, nCols};
}

function gridResolveFit(entry) {
  const model = document.getElementById("gridFitModelSelect").value;
  if (entry.fits && entry.fits[model]) return entry.fits[model];
  if (entry.model === model && entry.fit_x && entry.fit_y) {
    return {popt: entry.popt, fit_x: entry.fit_x, fit_y: entry.fit_y};
  }
  return null;
}

function gridCellTraces(process, quantity, methods, showFit) {
  const traces = [];
  const isMass = quantity.kind === "mass";
  methods.forEach(method => {
    const entry = getEntry(method, process, quantity);
    if (!entry) return;
    const color = gridColorForMethod(method);
    const methodLabel = DASHBOARD_DATA.methods[method].label;
    if (isMass) {
      const rec = (entry.definitions || {}).reco;
      if (rec && rec.edges && rec.edges.length) {
        traces.push(massHistTrace(rec, color, methodLabel));
      }
      if (showFit && entry.fit && entry.fit.fit_x && entry.fit.fit_y) {
        traces.push({
          x: entry.fit.fit_x, y: entry.fit.fit_y, mode: "lines", type: "scatter",
          name: methodLabel + " (fit)", line: {color: color, dash: "dash"},
          showlegend: false, hoverinfo: "skip",
        });
      }
    } else {
      const yKey = entry.sigma_over_E !== undefined ? "sigma_over_E" : "sigma";
      traces.push({
        x: entry.mid_points || [], y: entry[yKey] || [], mode: "markers", type: "scatter",
        name: methodLabel, marker: {color: color, size: 5},
      });
      const fit = showFit ? gridResolveFit(entry) : null;
      if (fit) {
        traces.push({
          x: fit.fit_x, y: fit.fit_y, mode: "lines", type: "scatter",
          name: methodLabel + " (fit)", line: {color: color, dash: "dash"},
          showlegend: false, hoverinfo: "skip",
        });
      }
    }
  });
  return traces;
}

function rebuildGridColorList() {
  const container = document.getElementById("gridColorList");
  container.innerHTML = "";
  gridSelectedMethods().forEach(method => {
    const row = document.createElement("div");
    row.className = "sel-row";
    const colorInput = document.createElement("input");
    colorInput.type = "color";
    colorInput.value = gridColorForMethod(method);
    colorInput.addEventListener("input", () => {
      gridMethodColors.set(method, colorInput.value);
      redrawGrid();
    });
    const label = DASHBOARD_DATA.methods[method].label;
    const span = document.createElement("span");
    span.textContent = label;
    span.title = label;
    row.appendChild(colorInput);
    row.appendChild(span);
    container.appendChild(row);
  });
}

function rebuildGridLegend(methods) {
  const container = document.getElementById("grid-legend");
  container.innerHTML = "";
  methods.forEach(method => {
    const item = document.createElement("div");
    item.className = "leg-item";
    const sw = document.createElement("div");
    sw.className = "leg-swatch";
    sw.style.background = gridColorForMethod(method);
    const span = document.createElement("span");
    span.textContent = DASHBOARD_DATA.methods[method].label;
    item.appendChild(sw);
    item.appendChild(span);
    container.appendChild(item);
  });
}

function redrawGrid() {
  rebuildGridColorList();
  const quantity = gridCurrentQuantity();
  const methods = gridSelectedMethods();
  const showFit = document.getElementById("gridShowFitCheckbox").checked;
  const {cellOf, nRows, nCols} = gridLayout();
  rebuildGridLegend(methods);

  const gridPlots = document.getElementById("gridPlots");
  gridPlots.style.gridTemplateColumns = `repeat(${Math.max(1, nCols)}, 1fr)`;
  gridPlots.innerHTML = "";

  const isMass = quantity.kind === "mass";
  const xAxisTitle = isMass ? "m_H [GeV]" :
    (quantity.kind === "eta_scan") ? "eta" :
    (quantity.kind === "costheta_scan") ? "cos(theta)" : "E_true [GeV]";
  const yAxisTitle = isMass ? "Events" :
    (quantity.kind === "angle") ? "sigma" : "sigma(E)/E";
  const xRange = isMass ? [90, 150] : null;

  for (let r = 0; r < nRows; r++) {
    for (let c = 0; c < nCols; c++) {
      const process = cellOf[r + "," + c];
      const cell = document.createElement("div");
      if (!process) {
        cell.className = "grid-cell empty";
        gridPlots.appendChild(cell);
        continue;
      }
      cell.className = "grid-cell";
      const title = document.createElement("div");
      title.className = "grid-cell-title";
      const procLabel = (DASHBOARD_DATA.process_meta[process] || {}).label || process;
      title.textContent = procLabel;
      title.title = procLabel;
      const plotDiv = document.createElement("div");
      plotDiv.className = "grid-cell-plot";
      cell.appendChild(title);
      cell.appendChild(plotDiv);
      gridPlots.appendChild(cell);

      const traces = gridCellTraces(process, quantity, methods, showFit);
      const layout = {
        margin: {t: 6, r: 6, b: 34, l: 46},
        showlegend: false,
        xaxis: {title: {text: xAxisTitle, font: {size: 9}}, tickfont: {size: 8}},
        yaxis: {title: {text: yAxisTitle, font: {size: 9}}, tickfont: {size: 8}},
        hovermode: "closest",
      };
      if (xRange) layout.xaxis.range = xRange;
      Plotly.react(plotDiv, traces, layout, {displayModeBar: false, responsive: true});
    }
  }
  saveStateToHash();
}

// --------------------------------------------------------------------------
// Grid presets: method comparisons mirroring the reference presentation.
//   - Detector / matching (slides 12 & 16): PF vs. Calo vs. PF+ideal matching.
//   - Higgs-mass clustering-algorithm scan (slide 19): Durham vs. anti-kt radii.
//   - Higgs-mass anti-kt energy recovery (slide 21): anti-kt radii with the
//     E-recovery correction (paired against the Durham baseline).
// --------------------------------------------------------------------------

function applyGridPreset(preset) {
  const methodSet = new Set(preset.methods);
  gridMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => {
    cb.checked = methodSet.has(cb.value);
  });
  preset.methods.forEach(m => {
    const color = preset.colorOf(m);
    if (color) gridMethodColors.set(m, color);
  });
  if (preset.quantity) gridQuantitySelect.value = preset.quantity;
  redrawGrid();
}

function buildGridPresets() {
  const presets = [];
  const durhamLike = methodNames.find(m => m === "PF_Durham") || methodNames[0];

  const detectorMethods = ["PF_Durham", "CaloJets_Durham", "PF_Durham_IdealMatching"]
    .filter(m => methodNames.includes(m));
  const detectorColors = {
    PF_Durham: "#1f77b4", CaloJets_Durham: "#2ca02c", PF_Durham_IdealMatching: "#ff7f0e",
  };

  const akMethods = methodNames
    .filter(m => antiKtRadius(m, false) !== null)
    .sort((a, b) => antiKtRadius(a, false) - antiKtRadius(b, false));
  const akRecMethods = methodNames
    .filter(m => antiKtRadius(m, true) !== null)
    .sort((a, b) => antiKtRadius(a, true) - antiKtRadius(b, true));

  const purpleFor = (list, m, radiusFn) => {
    const idx = list.indexOf(m);
    const lightness = 70 - (idx / Math.max(1, list.length - 1)) * 45;
    return hslToHex(270, 60, lightness);
  };

  // Slides 12 & 16: detector / matching comparison, jet energy resolution.
  if (detectorMethods.length > 1) {
    presets.push({
      title: "Detector / matching (energy resolution)",
      desc: "PF jets vs. Calo jets vs. PF + ideal matching, across all processes",
      methods: detectorMethods, colorOf: m => detectorColors[m], quantity: "energy:_all",
    });
    presets.push({
      title: "Detector / matching (Higgs mass)",
      desc: "PF jets vs. Calo jets vs. PF + ideal matching mH peaks, across all processes",
      methods: detectorMethods, colorOf: m => detectorColors[m], quantity: "mass",
    });
  }

  // Slide 19: Higgs mass, clustering-algorithm scan (Durham vs. anti-kt radii).
  if (durhamLike && akMethods.length) {
    const methods = [durhamLike, ...akMethods];
    presets.push({
      title: "Higgs mass: clustering-algorithm scan",
      desc: "Durham vs. anti-kt radii mH peaks, across all processes",
      methods,
      colorOf: m => (m === durhamLike ? "#1f77b4" : purpleFor(akMethods, m)),
      quantity: "mass",
    });
    presets.push({
      title: "Energy resolution: clustering-algorithm scan",
      desc: "Durham vs. anti-kt radii jet energy resolution, across all processes",
      methods,
      colorOf: m => (m === durhamLike ? "#1f77b4" : purpleFor(akMethods, m)),
      quantity: "energy:_all",
    });
  }

  // Slide 21: Higgs mass, anti-kt energy recovery (recovered radii vs. Durham).
  if (akRecMethods.length) {
    const methods = durhamLike ? [durhamLike, ...akRecMethods] : akRecMethods.slice();
    presets.push({
      title: "Higgs mass: anti-kt energy recovery",
      desc: "Durham vs. anti-kt radii with energy recovery, mH peaks across all processes",
      methods,
      // Matches slide 21: Durham in blue, AK-ER radii light->dark purple by radius.
      colorOf: m => (m === durhamLike ? "#1f77b4" : purpleFor(akRecMethods, m)),
      quantity: "mass",
    });
  }

  const container = document.getElementById("gridPresetList");
  container.innerHTML = "";
  presets.forEach(preset => {
    const btn = document.createElement("button");
    btn.className = "preset-btn";
    btn.innerHTML = `${preset.title}<span class="preset-desc">${preset.desc}</span>`;
    btn.addEventListener("click", () => applyGridPreset(preset));
    container.appendChild(btn);
  });
  if (presets.length === 0) {
    container.innerHTML = '<span style="font-size:11px;color:#999;">No presets available for this dataset.</span>';
  }
}

buildOptionList(gridMethodOptions, methodNames, n => DASHBOARD_DATA.methods[n].label, redrawGrid);
buildGridPresets();
gridQuantitySelect.addEventListener("change", redrawGrid);
document.getElementById("gridShowFitCheckbox").addEventListener("change", redrawGrid);
document.getElementById("gridFitModelSelect").addEventListener("change", redrawGrid);
document.getElementById("gridMethodsAllBtn").addEventListener("click", () => {
  gridMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = true; });
  redrawGrid();
});
document.getElementById("gridMethodsNoneBtn").addEventListener("click", () => {
  gridMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = false; });
  redrawGrid();
});
document.getElementById("gridMethodsAutoColorBtn").addEventListener("click", () => {
  const methods = gridSelectedMethods();
  gridMethodColors.clear();
  gridAutoColorIdx = 0;
  methods.forEach((m, i) => {
    const hue = (i * 360) / Math.max(1, methods.length);
    gridMethodColors.set(m, hslToHex(hue, 65, 45));
  });
  redrawGrid();
});

// --------------------------------------------------------------------------
// Tabs
// --------------------------------------------------------------------------

document.querySelectorAll(".tab-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".tab-content").forEach(c => c.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById("tab-" + btn.dataset.tab).classList.add("active");
    // Plotly can't size plots laid out while their container was display:none,
    // so (re)draw the grid the moment its tab becomes visible.
    if (btn.dataset.tab === "grid") redrawGrid();
    if (btn.dataset.tab === "coeff") redrawCoeff();
    saveStateToHash();
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
// Coefficients tab: how each fitted JER parameter (A, B, C, D...) trends across
// configurations - vs anti-kT radius, number of final-state jets, or clustering
// method. One panel per coefficient, error bars from the fit covariance. Fits
// live on the energy._all quantity (added by build_dashboard_data.py).
// --------------------------------------------------------------------------

// Coefficient display labels, in stored popt (function-signature) order, per
// base model. Kept explicit so the panels are labelled by physical meaning
// rather than a bare index. (The classic three/two-param popt is [stoch, const,
// noise] / [stoch, const].)
const COEFF_LABELS = {
  three_param: ["A: stochastic (1/&#8730;E)", "B: constant", "C: noise (1/E)"],
  three_param_quadrature: ["A: stochastic (1/&#8730;E)", "B: constant", "C: noise (1/E)"],
  two_param: ["A: stochastic (1/&#8730;E)", "B: constant"],
  two_param_quadrature: ["A: stochastic (1/&#8730;E)", "B: constant"],
  logE_4param: ["A: 1/&#8730;E", "B: constant", "C: 1/E", "D: log E"],
  invE15_noise: ["A: 1/&#8730;E", "B: constant", "C: 1/E", "E: 1/E^1.5"],
  invE2_noise: ["A: 1/&#8730;E", "B: constant", "C: 1/E", "E: 1/E^2"],
  logE_invE2_noStoch: ["B: constant", "C: 1/E", "E: 1/E^2", "D: log E"],
  logE_and_invE2: ["A: 1/&#8730;E", "B: constant", "C: 1/E", "E: 1/E^2", "D: log E"],
};
// Detailed explanations for each fit model, shown on the Coefficients tab for
// the selected method. Physics motivation + how it fared in the fit study
// (see src/plotting/fit_trials/; metrics are means over 159 method x process
// series, LOO CV = leave-one-out cross-validated MAE, the overfitting-robust
// figure of merit).
const FIT_MODEL_DESCRIPTIONS = {
  logE_4param_weighted: {
    title: "log(E) 4-parameter, error-weighted &mdash; the recommended fit",
    formula: "&sigma;/E = A/&#8730;E + B + C/E + D&middot;log(E)",
    body: [
      "The classic three-term resolution form (stochastic A/&#8730;E, constant B, noise C/E) plus a slow logarithmic drift D&middot;log(E). At jet level, calorimeter shower fluctuations, confusion/leakage in the clustering, and pile-up-like effects don't add up to a strictly constant high-energy floor; the log term absorbs that gentle residual energy dependence that a flat B cannot.",
      "It is fit weighting each point by its statistical uncertainty &delta;(&sigma;/E) = (&sigma;/E)/&#8730;(2N). In our study this was the <b>best-generalising model of all 45 tried</b>: its cross-validated error essentially equals its in-sample error, meaning it captures real structure rather than fitting noise. It beats the current production fit by ~3&times; on held-out points.",
    ],
    metrics: "LOO-CV MAE &asymp; 6.5&times;10&#8315;&#8308; &nbsp;|&nbsp; in-sample MAE &asymp; 6.8&times;10&#8315;&#8308; &nbsp;|&nbsp; rank #1 by CV",
  },
  logE_4param: {
    title: "log(E) 4-parameter (unweighted)",
    formula: "&sigma;/E = A/&#8730;E + B + C/E + D&middot;log(E)",
    body: [
      "Same functional form as the recommended fit, but with ordinary (unweighted) least squares, so every energy point counts equally in absolute terms. This gives a slightly better <em>in-sample</em> point-wise error than the weighted version, because unweighted least squares minimises exactly that quantity.",
      "Use this if you want the visually tightest curve through the plotted points and don't need the points weighted by their statistics. For a fit that generalises best to unseen energies, prefer the weighted variant.",
    ],
    metrics: "in-sample MAE &asymp; 4.5&times;10&#8315;&#8308; &nbsp;|&nbsp; LOO-CV MAE &asymp; 8.4&times;10&#8315;&#8308;",
  },
  invE2_noise: {
    title: "Four-term with an extra 1/E&sup2; (unweighted)",
    formula: "&sigma;/E = A/&#8730;E + B + C/E + e/E&sup2;",
    body: [
      "Adds a term even steeper than the usual 1/E noise term. The very-low-energy turn-up of the resolution is sharper than A/&#8730;E + C/E alone can follow; the e/E&sup2; term bends the curve up fast enough to track those first few GeV.",
      "Excellent in-sample accuracy (one of the best). Under cross-validation it is a touch behind the log(E) form, because the steep 1/E&sup2; term is more sensitive to the noisy lowest-energy point.",
    ],
    metrics: "in-sample MAE &asymp; 3.9&times;10&#8315;&#8308; &nbsp;|&nbsp; LOO-CV MAE &asymp; 8.2&times;10&#8315;&#8308;",
  },
  invE2_noise_weighted: {
    title: "Four-term with 1/E&sup2;, error-weighted",
    formula: "&sigma;/E = A/&#8730;E + B + C/E + e/E&sup2;",
    body: [
      "The 1/E&sup2; model fit with statistical weighting. Weighting pulls the fit toward the high-statistics mid-energy bins and slightly tames the sensitivity of the steep low-E term.",
      "A solid, physically-motivated alternative to the log(E) fit when you specifically want a steep low-energy term rather than a high-energy drift term.",
    ],
    metrics: "in-sample MAE &asymp; 4.5&times;10&#8315;&#8308; &nbsp;|&nbsp; LOO-CV MAE &asymp; 8.8&times;10&#8315;&#8308;",
  },
  invE15_noise: {
    title: "Four-term with a milder 1/E^1.5",
    formula: "&sigma;/E = A/&#8730;E + B + C/E + e/E^1.5",
    body: [
      "A gentler cousin of the 1/E&sup2; model: the extra low-energy term falls as 1/E^1.5, between the noise term (1/E) and the steep 1/E&sup2;. It captures the turn-up while being less dominated by the single lowest-energy point.",
      "Among the best on both in-sample and cross-validated error &mdash; a good compromise if the 1/E&sup2; term looks too aggressive for your data.",
    ],
    metrics: "in-sample MAE &asymp; 4.0&times;10&#8315;&#8308; &nbsp;|&nbsp; LOO-CV MAE &asymp; 8.1&times;10&#8315;&#8308;",
  },
  logE_invE2_noStoch: {
    title: "No stochastic term: constant + 1/E + 1/E&sup2; + log(E)",
    formula: "&sigma;/E = B + C/E + e/E&sup2; + D&middot;log(E)",
    body: [
      "Drops the textbook A/&#8730;E stochastic term entirely, relying on the steep (1/E, 1/E&sup2;) and slow (log E) terms instead. Included as a deliberate test of whether the &#8730;E term actually earns its place once the other terms are present.",
      "Remarkably, it fits nearly as well as the models that keep the &#8730;E term &mdash; a hint that at jet level (as opposed to single particles) the classic 1/&#8730;E scaling is not strongly preferred by these data. Useful as a physics cross-check rather than a default.",
    ],
    metrics: "in-sample MAE &asymp; 4.1&times;10&#8315;&#8308; &nbsp;|&nbsp; LOO-CV MAE &asymp; 8.3&times;10&#8315;&#8308;",
  },
  logE_and_invE2: {
    title: "Five-parameter: log(E) + 1/E&sup2; (unweighted)",
    formula: "&sigma;/E = A/&#8730;E + B + C/E + e/E&sup2; + D&middot;log(E)",
    body: [
      "Combines both winning extra terms (the steep 1/E&sup2; and the slow log E) on top of the classic three. With five free parameters it achieves the <b>best in-sample accuracy of any model tried</b>.",
      "<b>Caveat:</b> it <em>overfits</em>. Under leave-one-out cross-validation it drops well down the ranking &mdash; its held-out error is ~3&times; worse than its in-sample error, because the 5th parameter mostly fits noise on the 9&ndash;13-point series. Shown so you can see the trade-off; not recommended as the default.",
    ],
    metrics: "in-sample MAE &asymp; 2.5&times;10&#8315;&#8308; (best) &nbsp;|&nbsp; LOO-CV MAE &asymp; 8.9&times;10&#8315;&#8308; (overfits)",
  },
  logE_and_invE2_weighted: {
    title: "Five-parameter log(E) + 1/E&sup2;, error-weighted",
    formula: "&sigma;/E = A/&#8730;E + B + C/E + e/E&sup2; + D&middot;log(E)",
    body: [
      "The 5-parameter model with statistical weighting. Like its unweighted twin it has superb in-sample accuracy but overfits under cross-validation.",
      "Included for completeness / comparison. Prefer a 4-parameter model unless you have a specific reason to want all five terms.",
    ],
    metrics: "in-sample MAE &asymp; 2.7&times;10&#8315;&#8308; &nbsp;|&nbsp; LOO-CV MAE &asymp; 9.6&times;10&#8315;&#8308; (overfits)",
  },
  three_param: {
    title: "Classic 3-parameter (unbounded)",
    formula: "&sigma;/E = A/&#8730;E + B + C/E",
    body: [
      "The textbook calorimeter resolution parameterisation: stochastic (A/&#8730;E), constant (B) and noise (C/E) terms added linearly. This is the historical default for jet energy resolution.",
      "Fit here with no parameter bounds. Even this small change (the production version bounds B to [0.005, 0.04]) already roughly halves the error versus the bounded incumbent &mdash; but it is still clearly beaten by the four-term models above, which have the extra freedom to follow the low-E turn-up or high-E drift.",
    ],
    metrics: "in-sample MAE &asymp; 8.4&times;10&#8315;&#8308; &nbsp;|&nbsp; LOO-CV MAE &asymp; 1.2&times;10&#8315;&#179;",
  },
  three_param_weighted: {
    title: "Classic 3-parameter, error-weighted",
    formula: "&sigma;/E = A/&#8730;E + B + C/E",
    body: [
      "The classic three-term form fit with statistical weighting (Option A) and the production B-bound. Weighting makes the (statistically meaningful) &chi;&sup2; smaller, but with only three terms the model still cannot follow the low-energy turn-up as well as the four-term fits.",
      "This is the closest apples-to-apples 'improved incumbent' &mdash; useful as the reference the newer forms should be compared against.",
    ],
    metrics: "3-term family; clearly improved on by the 4-term models",
  },
  two_param: {
    title: "2-parameter (stochastic + constant)",
    formula: "&sigma;/E = A/&#8730;E + B",
    body: [
      "The simplest form: just a stochastic term and a constant floor, with no noise term. Appropriate when the low-energy noise contribution is negligible.",
      "Fewest parameters, so the most robust against overfitting, but it systematically misses the low-energy turn-up here. Best used as a sanity-check baseline.",
    ],
    metrics: "in-sample MAE &asymp; 2.2&times;10&#8315;&#179;",
  },
  two_param_weighted: {
    title: "2-parameter, error-weighted",
    formula: "&sigma;/E = A/&#8730;E + B",
    body: [
      "The two-parameter form with statistical weighting. Same limitations as the unweighted version &mdash; too few terms to follow the full energy dependence &mdash; but a useful minimal reference.",
    ],
    metrics: "2-term baseline",
  },
  staged_const_first: {
    title: "3-parameter staged: constant pinned to the high-E tail",
    formula: "&sigma;/E = A/&#8730;E + B + C/E, &nbsp; B = mean(&sigma;/E over the 3 highest-E points)",
    body: [
      "The classic three-term form, but fitted in two stages rather than all at once. At the top of the energy range the stochastic and noise terms have died away, so &sigma;/E there is essentially the constant floor: stage 1 reads B off as the mean of the three highest-energy points, stage 2 holds B fixed and least-squares fits only A and C to every point. The appeal is that B gets a direct physical reading instead of being traded off against the other two terms by the fitter.",
      "In practice it is the <b>worst</b> of the strategies tried: the highest-energy points are the lowest-statistics ones, so B inherits their noise and, being frozen, that error propagates into A and C across the whole range. Releasing B afterwards (using the tail only as a starting guess) recovers the ordinary unbounded fit exactly, which is the better way to use the idea. Shown here so the effect is visible per method and process rather than only as a leaderboard row.",
    ],
    metrics: "in-sample MAE &asymp; 2.7&times;10&#8315;&#179; &nbsp;|&nbsp; LOO-CV MAE &asymp; 3.0&times;10&#8315;&#179; &nbsp;|&nbsp; last of 45 by CV",
  },
};

function renderCoeffModelDesc(model) {
  const d = FIT_MODEL_DESCRIPTIONS[model];
  const el = document.getElementById("coeffModelDesc");
  if (!d) { el.innerHTML = ""; return; }
  el.innerHTML = `<h3>${d.title}</h3>`
    + `<p><span class="formula">${d.formula}</span></p>`
    + d.body.map(p => `<p>${p}</p>`).join("")
    + (d.metrics ? `<div class="metrics"><b>Fit study:</b> ${d.metrics}</div>` : "");
}

const COEFF_ENERGY_Q = {kind: "energy", part: "_all"};
const COEFF_X_HINT = {
  radius: "anti-kT radii only (Durham/Calo have no radius and are skipped). One line per process; dotted = with energy recovery.",
  njets: "Points are processes grouped by their final-state jet count; one colour per method.",
  method: "One point per clustering method; one colour per process.",
};
// stable per-method colours for the njets view (series = method)
const coeffMethodColor = {};
methodNames.forEach((m, i) => { coeffMethodColor[m] = AUTO_COLORS[i % AUTO_COLORS.length]; });

function coeffLabelsFor(model) { return COEFF_LABELS[fitBaseModelOf(model)] || null; }
function coeffRadius(method) { const m = method.match(/AntiKtR(\d+)/); return m ? parseInt(m[1], 10) / 10 : null; }
function coeffIsRecovery(method) { return /E_recovery/i.test(method); }
function methodLabelOf(m) { return (DASHBOARD_DATA.methods[m] || {}).label || m; }
function procLabelOf(p) { return (DASHBOARD_DATA.process_meta[p] || {}).label || p; }
function procColorOf(p) { return (DASHBOARD_DATA.process_meta[p] || {}).color || AUTO_COLORS[0]; }

const coeffFitModelSelect = document.getElementById("coeffFitModelSelect");
const coeffXaxisSelect = document.getElementById("coeffXaxisSelect");
const coeffMethodOptions = document.getElementById("coeffMethodOptions");
const coeffProcessOptions = document.getElementById("coeffProcessOptions");
populateFitModelSelect(coeffFitModelSelect, "logE_4param_weighted");
buildOptionList(coeffMethodOptions, methodNames, methodLabelOf, redrawCoeff);
buildOptionList(coeffProcessOptions, processNames, procLabelOf, redrawCoeff);
// start with everything selected so the tab shows something immediately
coeffMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = true; });
coeffProcessOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = true; });

// Compact tick label for the (categorical) clustering-method x-axis, so 15
// long method names don't force the panel wider than its container.
function shortMethodLabel(m) {
  if (m === "CaloJets_Durham") return "Calo";
  if (m === "PF_Durham") return "PF Durham";
  if (m === "PF_Durham_IdealMatching") return "PF Ideal";
  const rr = m.match(/AntiKtR(\d+)/);
  if (rr) return "AK" + rr[1] + (coeffIsRecovery(m) ? " ER" : "");
  return m;
}
// Sort key so the method axis reads left-to-right sensibly: radius-less methods
// (Durham/Calo) first, then anti-kT by increasing radius, E-recovery just after
// its plain counterpart.
function coeffMethodRank(m) {
  const r = coeffRadius(m);
  if (r === null) return -1;
  return r + (coeffIsRecovery(m) ? 0.001 : 0);
}

// Map one point to its series (grouping/colour/x) for the chosen x-axis.
function coeffSeriesOf(p, xaxis) {
  if (xaxis === "radius") {
    return {key: p.process + (p.rec ? "|rec" : "|std"),
            name: procLabelOf(p.process) + (p.rec ? " (E-rec)" : ""),
            color: procColorOf(p.process), dash: p.rec ? "dot" : "solid", x: p.radius};
  }
  if (xaxis === "njets") {
    return {key: p.method, name: methodLabelOf(p.method),
            color: coeffMethodColor[p.method] || AUTO_COLORS[0], dash: "solid", x: p.njets};
  }
  return {key: p.process, name: procLabelOf(p.process),
          color: procColorOf(p.process), dash: "solid", x: shortMethodLabel(p.method)};
}

function redrawCoeff() {
  const model = coeffFitModelSelect.value;
  const xaxis = coeffXaxisSelect.value;
  const methods = getChecked(coeffMethodOptions);
  const processes = getChecked(coeffProcessOptions);
  const labels = coeffLabelsFor(model);
  const plotsDiv = document.getElementById("coeffPlots");
  const legendDiv = document.getElementById("coeffLegend");
  const msg = document.getElementById("coeffMsg");
  document.getElementById("coeffXaxisHint").textContent = COEFF_X_HINT[xaxis] || "";
  renderCoeffModelDesc(model);
  plotsDiv.innerHTML = ""; legendDiv.innerHTML = ""; msg.textContent = "";

  const pts = [];
  methods.forEach(method => processes.forEach(process => {
    const entry = getEntry(method, process, COEFF_ENERGY_Q);
    const fit = entry && entry.fits && entry.fits[model];
    if (!fit || !fit.popt) return;
    pts.push({method, process, popt: fit.popt, err: fit.popt_err || [],
              radius: coeffRadius(method), rec: coeffIsRecovery(method),
              njets: (DASHBOARD_DATA.process_meta[process] || {}).n_jets});
  }));
  if (!pts.length) { msg.textContent = "No fits available for this fit model / selection."; return; }

  const nCoef = labels ? labels.length : Math.max.apply(null, pts.map(p => p.popt.length));
  const legendMap = new Map();  // name -> {color, dash}
  let anyTrace = false;

  // Per-mode x-axis config. `automargin` lets Plotly grow the margin (not the
  // width) for tick labels, and a fixed width-fitting layout keeps every panel
  // inside its grid cell so there's never a horizontal scrollbar. Rebuilt on
  // every redraw, so the scale re-fits after any dropdown/checkbox change.
  let xAxisLayout;
  if (xaxis === "method") {
    const catArray = methods.slice().sort((a, b) => coeffMethodRank(a) - coeffMethodRank(b))
      .map(shortMethodLabel);
    xAxisLayout = {title: {text: "", font: {size: 11}}, type: "category",
                   categoryorder: "array", categoryarray: catArray,
                   tickangle: -45, automargin: true, tickfont: {size: 10}};
  } else if (xaxis === "njets") {
    xAxisLayout = {title: {text: "# final-state jets", font: {size: 11}}, type: "linear",
                   tickvals: [2, 4, 6], range: [1, 7], automargin: true};
  } else {  // radius
    xAxisLayout = {title: {text: "anti-kT radius R", font: {size: 11}}, type: "linear",
                   dtick: 0.2, autorange: true, automargin: true};
  }

  for (let ci = 0; ci < nCoef; ci++) {
    const groups = new Map();
    pts.forEach(p => {
      if (p.popt[ci] === undefined) return;
      const s = coeffSeriesOf(p, xaxis);
      if (s.x === null || s.x === undefined) return;
      if (!groups.has(s.key)) groups.set(s.key, {s, xs: [], ys: [], es: []});
      const g = groups.get(s.key);
      g.xs.push(s.x); g.ys.push(p.popt[ci]); g.es.push((p.err && p.err[ci]) || 0);
      legendMap.set(s.name, {color: s.color, dash: s.dash});
    });
    const traces = [];
    groups.forEach(g => {
      let order = g.xs.map((_, i) => i);
      if (typeof g.xs[0] === "number") order.sort((a, b) => g.xs[a] - g.xs[b]);
      traces.push({
        x: order.map(i => g.xs[i]), y: order.map(i => g.ys[i]),
        mode: xaxis === "radius" ? "lines+markers" : "markers", type: "scatter",
        name: g.s.name, showlegend: false,
        marker: {color: g.s.color, size: 7}, line: {color: g.s.color, dash: g.s.dash},
        error_y: {type: "data", array: order.map(i => g.es[i]), visible: true,
                  thickness: 1, width: 3, color: g.s.color},
      });
    });
    if (traces.length) anyTrace = true;
    const cell = document.createElement("div");
    cell.style.cssText = "border:1px solid #eee; border-radius:4px; background:#fff; padding:2px; min-width:0; overflow:hidden;";
    plotsDiv.appendChild(cell);
    Plotly.react(cell, traces, {
      autosize: true,
      margin: {t: 30, l: 54, r: 12, b: 44},
      title: {text: labels ? labels[ci] : "param " + ci, font: {size: 12}},
      xaxis: JSON.parse(JSON.stringify(xAxisLayout)),
      yaxis: {title: {text: "coefficient", font: {size: 10}}, zeroline: true, automargin: true},
      height: 300, showlegend: false, hovermode: "closest",
    }, {responsive: true, displayModeBar: false});
    // Container just changed size/visibility; force Plotly to re-fit to it so a
    // panel can never end up wider than its cell (the cause of the h-scroll).
    Plotly.Plots.resize(cell);
  }

  if (!anyTrace) {
    msg.textContent = xaxis === "radius"
      ? "None of the selected methods have an anti-kT radius (Durham/Calo are radius-less)."
      : "Nothing to plot for this selection.";
  }
  // shared legend
  legendMap.forEach((v, name) => {
    const item = document.createElement("span");
    item.style.cssText = "display:flex; align-items:center; gap:5px;";
    const sw = document.createElement("span");
    sw.style.cssText = `width:16px; height:3px; border-radius:2px; background:${v.color};`
      + (v.dash === "dot" ? " outline:1px dotted " + v.color + "; background:transparent;" : "");
    item.appendChild(sw);
    item.appendChild(document.createTextNode(name));
    legendDiv.appendChild(item);
  });
}

coeffFitModelSelect.addEventListener("change", redrawCoeff);
coeffXaxisSelect.addEventListener("change", redrawCoeff);
document.getElementById("coeffMethodsAllBtn").addEventListener("click", () => {
  coeffMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = true; });
  redrawCoeff();
});
document.getElementById("coeffMethodsNoneBtn").addEventListener("click", () => {
  coeffMethodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = false; });
  redrawCoeff();
});
document.getElementById("coeffProcessesAllBtn").addEventListener("click", () => {
  coeffProcessOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = true; });
  redrawCoeff();
});
document.getElementById("coeffProcessesNoneBtn").addEventListener("click", () => {
  coeffProcessOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = false; });
  redrawCoeff();
});
function coeffSelectByNJets(n) {
  coeffProcessOptions.querySelectorAll("input[type=checkbox]").forEach(cb => {
    cb.checked = ((DASHBOARD_DATA.process_meta[cb.value] || {}).n_jets === n);
  });
  redrawCoeff();
}
document.getElementById("coeffProcesses2jBtn").addEventListener("click", () => coeffSelectByNJets(2));
document.getElementById("coeffProcesses4jBtn").addEventListener("click", () => coeffSelectByNJets(4));
document.getElementById("coeffProcesses6jBtn").addEventListener("click", () => coeffSelectByNJets(6));

// Keep the coefficient panels fitted to their container when the window is
// resized while the tab is open (Plotly's responsive handler only refits plots
// that were visible at creation; these are (re)built on demand).
window.addEventListener("resize", () => {
  if (!document.getElementById("tab-coeff").classList.contains("active")) return;
  document.querySelectorAll("#coeffPlots > div").forEach(c => Plotly.Plots.resize(c));
});

// --------------------------------------------------------------------------
// Statistics tab: fit coefficients, event counts, filter pass rates.
// --------------------------------------------------------------------------

// One "value +/- error" term, error omitted when unavailable (e.g. angular
// fits that weren't reweighted). popt/popt_err share the raw curve_fit order
// [A, C, B] for the three-param model, hence the index gymnastics below.
function pmTerm(name, val, err) {
  const v = Number(val).toFixed(3);
  if (err === undefined || err === null || !isFinite(err)) return `${name}=${v}`;
  return `${name}=${v}±${Number(err).toFixed(3)}`;
}

// Param labels in popt (function-signature) order for the experiment-derived
// energy-JER models added in build_dashboard_data.py. The classic three/two
// param fits keep their historical A/B/C display order handled separately below.
const FIT_PARAM_NAMES = {
  logE_4param: ["A", "B", "C", "D_log"],
  invE15_noise: ["A", "B", "C", "E_1.5"],
  invE2_noise: ["A", "B", "C", "E_2"],
  logE_invE2_noStoch: ["B", "C", "E_2", "D_log"],
  logE_and_invE2: ["A", "B", "C", "E_2", "D_log"],
};

function formatFitParams(model, popt, popt_err) {
  const e = popt_err || [];
  const base = fitBaseModelOf(model);
  if (base === "two_param" || base === "two_param_quadrature") {
    return `${pmTerm("A", popt[0], e[0])} ${pmTerm("C", popt[1], e[1])}`;
  }
  if (base === "three_param" || base === "three_param_quadrature") {
    // Classic form's stored popt order is [A, C, B]; display A, B, C.
    return `${pmTerm("A", popt[0], e[0])} ${pmTerm("B", popt[2], e[2])} ${pmTerm("C", popt[1], e[1])}`;
  }
  if (FIT_PARAM_NAMES[base]) {
    // Experiment-derived models: popt is already in signature/label order.
    return FIT_PARAM_NAMES[base].map((nm, i) => pmTerm(nm, popt[i], e[i])).join(" ");
  }
  if (model === "gaussian") {
    return `mu=${popt[0].toFixed(2)} sigma=${popt[1].toFixed(2)}`;
  }
  return popt.map(v => Number(v).toFixed(3)).join(", ");
}

function formatChi2(fit) {
  return (fit && fit.chi2_ndf !== undefined && fit.chi2_ndf !== null)
    ? Number(fit.chi2_ndf).toFixed(2) : "";
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
      const title = isObj && cell.title ? ` title="${cell.title}"` : "";
      html += `<td${cls}${title}>${text === undefined || text === null ? "" : text}</td>`;
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
            rows.push([methodLabel, procLabel, "Higgs mass peak", fit.model, formatFitParams(fit.model, fit.popt), ""]);
          }
          return;
        }
        const entry = getEntry(method, process, q);
        if (!entry) return;
        if (entry.fits) {
          Object.keys(entry.fits).forEach(model => {
            const fit = entry.fits[model];
            if (fit) rows.push([methodLabel, procLabel, q.label, model,
              formatFitParams(model, fit.popt, fit.popt_err), formatChi2(fit)]);
          });
        } else if (entry.model && entry.popt) {
          rows.push([methodLabel, procLabel, q.label, entry.model,
            formatFitParams(entry.model, entry.popt, entry.popt_err), ""]);
        }
      });
    });
  });
  buildStatsTable(document.getElementById("fitCoeffTable"),
    ["Method", "Process", "Quantity", "Model", "Parameters", "χ²/ndf"], rows);
}

function renderCountAndPassRateTables() {
  const stats = DASHBOARD_DATA.stats || {};
  const processes = Object.keys(stats).sort();
  const folderSet = new Set();
  processes.forEach(p => Object.keys(stats[p]).forEach(f => folderSet.add(f)));
  const folders = Array.from(folderSet).sort();

  const compactFmt = new Intl.NumberFormat(undefined, {notation: "compact", maximumFractionDigits: 1});
  const countRows = processes.map(p => {
    const row = [p];
    folders.forEach(f => {
      const m = stats[p][f];
      if (!m) { row.push(""); return; }
      row.push({
        text: `${compactFmt.format(m.before)} / ${compactFmt.format(m.after)}`,
        title: `${Math.round(m.before).toLocaleString()} / ${Math.round(m.after).toLocaleString()}`,
      });
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
document.getElementById("showFitCheckbox").addEventListener("change", redraw);
document.getElementById("showErrorBarsCheckbox").addEventListener("change", redraw);
document.getElementById("fitModelSelect").addEventListener("change", redraw);

document.getElementById("showFullHistCheckbox").addEventListener("change", () => {
  redraw();
  redrawHistPlot();
});

document.getElementById("clearHistBtn").addEventListener("click", () => {
  selectedHistPoints.clear();
  redrawHistPlot();
});
document.getElementById("normalizeHistCheckbox").addEventListener("change", redrawHistPlot);

// --------------------------------------------------------------------------
// Bulk-selection convenience buttons (Methods/Processes All/None, jet-count
// quick-select) and "Auto colors" to drop back to the default coloring
// (process_config.py colors / auto-cycled fallback) for every combo.
// --------------------------------------------------------------------------

document.getElementById("methodsAllBtn").addEventListener("click", () => {
  methodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = true; });
  redraw();
});
document.getElementById("methodsNoneBtn").addEventListener("click", () => {
  methodOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = false; });
  redraw();
});
document.getElementById("processesAllBtn").addEventListener("click", () => {
  processOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = true; });
  redraw();
});
document.getElementById("processesNoneBtn").addEventListener("click", () => {
  processOptions.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = false; });
  redraw();
});

function selectProcessesByNJets(n) {
  processOptions.querySelectorAll("input[type=checkbox]").forEach(cb => {
    const meta = DASHBOARD_DATA.process_meta[cb.value] || {};
    cb.checked = meta.n_jets === n;
  });
  redraw();
}
document.getElementById("processes2jBtn").addEventListener("click", () => selectProcessesByNJets(2));
document.getElementById("processes4jBtn").addEventListener("click", () => selectProcessesByNJets(4));
document.getElementById("processes6jBtn").addEventListener("click", () => selectProcessesByNJets(6));

document.getElementById("autoColorsBtn").addEventListener("click", () => {
  // Assign evenly-spaced hues across the currently selected combos, so they
  // are guaranteed visually distinct from each other - unlike falling back to
  // process_config.py colors, which can make two different methods showing
  // the same process look identical.
  const methods = getChecked(methodOptions);
  const processes = getChecked(processOptions);
  const combos = [];
  methods.forEach(m => processes.forEach(p => combos.push(m + "||" + p)));

  selectedColors.clear();
  autoColorIdx = 0;
  combos.forEach((k, i) => {
    const hue = (i * 360) / Math.max(1, combos.length);
    selectedColors.set(k, hslToHex(hue, 65, 45));
  });
  redraw();
});

// Restore state from the URL hash if present (bookmarked/shared link);
// otherwise fall back to sensible defaults so the dashboard isn't empty.
const initialState = decodeState(location.hash.replace(/^#/, ""));
if (!initialState) {
  if (methodNames.length) methodOptions.querySelector("input[type=checkbox]").checked = true;
  if (processNames.length) processOptions.querySelector("input[type=checkbox]").checked = true;
}

// Restore state on browser back/forward navigation too. history.replaceState
// (used by saveStateToHash) never fires "hashchange", so this only reacts to
// actual user navigation, not our own writes.
window.addEventListener("hashchange", () => {
  const state = decodeState(location.hash.replace(/^#/, ""));
  if (state) applyState(state);
});

Plotly.newPlot("mainPlot", [], {margin: {t: 20}}).then(() => {
  attachClickHandler();
  if (initialState) {
    applyState(initialState);
  } else {
    redraw();
  }
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
