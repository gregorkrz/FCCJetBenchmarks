"""Stage 4 of the resolution pipeline (no ROOT, no fitting — pure packaging).

Reads the consolidated `dashboard_data.json` produced by
`build_dashboard_data.py` and writes a single self-contained HTML file
(Plotly.js loaded from CDN, vanilla JS, no build step) that lets the user
interactively explore the resolution grid:

- Multi-select filters for method (jet clustering algo / detector config),
  process, and quantity (jet-part energy resolution, angular resolution,
  eta/cos(theta) scans).
- A color picker per selected (method, process) combination, defaulting to
  the colors from `process_config.py` with automatic fallback colors.
- A main Plotly scatter+line plot of resolution vs. energy (or eta/costheta)
  for every selected combination, including the fitted curve.
- Clicking on a data point shows the underlying per-bin histogram (and the
  low/high/MPV fit values as vertical lines) in a secondary plot.

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
  #histHint { font-size: 12px; color: #777; padding: 4px 0; }
</style>
</head>
<body>
<h1>FCC Jet Resolution Dashboard</h1>
<div id="layout">
  <div id="controls">
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
    <div id="histHint">Click on a data point above to show the histogram used to compute that resolution value.</div>
    <div id="histPlot"></div>
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

function redraw() {
  const methods = getChecked(methodOptions);
  const processes = getChecked(processOptions);
  const quantity = currentQuantity();

  const combos = [];
  methods.forEach(m => processes.forEach(p => combos.push({method: m, process: p})));
  rebuildSelectionList(combos);

  const traces = [];
  currentPointIndex = [];

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

function attachClickHandler() {
  const gd = document.getElementById("mainPlot");
  gd.on("plotly_click", function(data) {
    const pt = data.points[0];
    const traceIdx = pt.curveNumber;
    const info = currentPointIndex[traceIdx];
    if (!info) return;
    const binIdx = pt.pointIndex;
    showHistogram(info, binIdx);
  });
}

function showHistogram(info, binIdx) {
  const bins = (info.entry.bins || []);
  const rec = bins[binIdx];
  if (!rec) {
    Plotly.purge("histPlot");
    return;
  }
  const edges = rec.edges || [];
  const y = rec.y || [];
  const centers = [];
  for (let i = 0; i < edges.length - 1; i++) {
    centers.push(0.5 * (edges[i] + edges[i + 1]));
  }
  const widths = [];
  for (let i = 0; i < edges.length - 1; i++) {
    widths.push(edges[i + 1] - edges[i]);
  }

  const histTrace = {
    x: centers,
    y: y,
    type: "bar",
    width: widths,
    marker: {color: "#888"},
    name: "histogram",
  };

  const yMax = Math.max(1, ...y);
  const vline = (xval, name, dash) => ({
    x: [xval, xval], y: [0, yMax], mode: "lines", type: "scatter",
    name: name, line: {color: "red", dash: dash || "solid"},
  });

  const traces = [histTrace];
  if (rec.low !== undefined) traces.push(vline(rec.low, "low", "dot"));
  if (rec.high !== undefined) traces.push(vline(rec.high, "high", "dot"));
  if (rec.mpv !== undefined) traces.push(vline(rec.mpv, "MPV", "solid"));

  const title = `${DASHBOARD_DATA.methods[info.method].label} / ` +
    `${(DASHBOARD_DATA.process_meta[info.process] || {}).label || info.process} ` +
    `[${rec.lo}, ${rec.hi}] GeV (N=${rec.n_jets || "?"})`;

  Plotly.react("histPlot", traces, {
    title: title,
    margin: {t: 40},
    showlegend: true,
  }, {responsive: true});
}

buildOptionList(methodOptions, methodNames, n => DASHBOARD_DATA.methods[n].label, redraw);
buildOptionList(processOptions, processNames, n => (DASHBOARD_DATA.process_meta[n] || {}).label || n, redraw);
quantitySelect.addEventListener("change", redraw);

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
