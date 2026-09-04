"""Regenerate the presentation's JER figures with a chosen energy-dependence fit.

The slides ("Jet Energy Resolution (Durham PF)", the 5x3 "full plots" grids, the
per-flavour overlays) all show the production `three_param` fit
`S/sqrt(E) + C + N/E`. This script redraws the same figures with any fit model
stored in `dashboard_data.json`; the default is `logE_4param_weighted`
(`S/sqrt(E) + C + N/E + D*log(E)`, stat-error weighted) - the best-generalising
form found by the `fit_trials/` study (leave-one-out CV MAE 6.5e-4 vs 2.0e-3 for
the production fit; see fit_trials/README.md).

No ROOT access or refitting is needed: `build_dashboard_data.py` already stores
every candidate fit (popt, curve, popt_err, chi2/ndf) per method x process, so
this only reads the JSON and draws.

Also emits the two Higgs-mass overlays (b-jet / light-flavour final states) as
single-panel figures, since they belong in the same slide set - the fit model
does not affect those.

Usage:
    source env.sh
    python3 src/plotting/presentation_jer_plots.py                  # logE_4param_weighted
    python3 src/plotting/presentation_jer_plots.py --model three_param
    python3 src/plotting/presentation_jer_plots.py --list-models
"""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Config mirroring the presentation figures
# ---------------------------------------------------------------------------

DEFAULT_MODEL = "logE_4param_weighted"

METHOD_LABELS = {
    "PF_Durham": "PFJets",
    "PF_Durham_IdealMatching": "PFJets + Ideal Matching",
    "CaloJets_Durham": "CaloJets",
}
METHOD_COLORS = {
    "PF_Durham": "blue",
    "PF_Durham_IdealMatching": "orange",
    "CaloJets_Durham": "green",
}

# The generalized e+e- anti-kt radius scan, with and without energy recovery.
AK_RADII = ["04", "06", "08", "10", "12", "14"]
AK_METHODS = [f"PF_AntiKtR{r}" for r in AK_RADII]
AK_ER_METHODS = [f"PF_E_recovery_AntiKtR{r}" for r in AK_RADII]
for _i, _r in enumerate(AK_RADII):
    _c = plt.cm.plasma(0.05 + 0.8 * _i / (len(AK_RADII) - 1))
    METHOD_LABELS[f"PF_AntiKtR{_r}"] = f"ee-AK R={int(_r) / 10:.1f}"
    METHOD_COLORS[f"PF_AntiKtR{_r}"] = _c
    METHOD_LABELS[f"PF_E_recovery_AntiKtR{_r}"] = f"ee-AK R={int(_r) / 10:.1f} (E-rec)"
    METHOD_COLORS[f"PF_E_recovery_AntiKtR{_r}"] = _c

# The two 2-jet processes shown on the "Jet Energy Resolution (Durham PF)" slide.
SELECTED_PROCESSES = ["p8_ee_ZH_vvqq_ecm240", "p8_ee_ZH_vvgg_ecm240"]

# Which slide each output PDF corresponds to (written into the output README).
FIGURE_DESCRIPTIONS = {
    "JER_grid_PFJets_vs_IdealMatching.pdf":
        'slide "Jet Energy Resolution (full plots)" - PFJets vs PFJets + Ideal Matching',
    "JER_grid_PFJets_CaloJets_IdealMatching.pdf":
        'slide "Jet energy resolution: comparison with calo jets (full plots)"',
    "JER_selected_2jet_processes.pdf":
        'slide "Jet Energy Resolution (Durham PF)" - the two Z(->vv)H panels',
    "JER_overlay_bjet_processes.pdf":
        'slide "Jet Energy Resolution (Durham PF)" - b-jet final states + response',
    "JER_overlay_lightflavour_processes.pdf":
        "same overlay for the light-flavour final states",
    "JER_grid_fit_comparison_PFJets.pdf":
        "new: both fit models overlaid per process (PFJets), with chi2/ndf",
    "JER_grid_fit_comparison_CaloJets.pdf":
        'new: same for CaloJets (the slide with the "issues with fitting" caveat)',
    "JER_fit_residuals.pdf":
        "new: fit - data residuals vs E for both models, all processes overlaid",
    "mH_bjet_processes.pdf":
        'slide "Results: Durham algorithm with PFO" - mH, b-jet final states (single panel)',
    "mH_lightflavour_processes.pdf":
        "same for the light-flavour final states (single panel)",
}

# Display labels for the fitted parameters, keyed by base model. Order follows
# the function signature; the printed order is S, N, C, extras (as on the
# slides, where the noise term N is printed before the constant C).
PARAM_DISPLAY = {
    "three_param": ["S", "C", "N"],
    "three_param_quadrature": ["S", "C", "N"],
    "two_param": ["S", "C"],
    "two_param_quadrature": ["S", "C"],
    "logE_4param": ["S", "C", "N", "D"],
    "invE15_noise": ["S", "C", "N", "E15"],
    "invE2_noise": ["S", "C", "N", "E2"],
    "logE_invE2_noStoch": ["C", "N", "E2", "D"],
    "logE_and_invE2": ["S", "C", "N", "E2", "D"],
}
PRINT_ORDER = ["S", "N", "C", "D", "E15", "E2"]
# Small-magnitude coefficients need more digits than the slides' round(x, 2).
SCIENTIFIC = {"D", "E15", "E2"}

MODEL_FORMULA = {
    "three_param": r"$S/\sqrt{E} + C + N/E$",
    "three_param_quadrature": r"$S/\sqrt{E} \oplus C \oplus N/E$",
    "logE_4param": r"$S/\sqrt{E} + C + N/E + D\ln E$",
    "invE15_noise": r"$S/\sqrt{E} + C + N/E + E_{1.5}/E^{1.5}$",
    "invE2_noise": r"$S/\sqrt{E} + C + N/E + E_{2}/E^{2}$",
    "logE_invE2_noStoch": r"$C + N/E + E_{2}/E^{2} + D\ln E$",
    "logE_and_invE2": r"$S/\sqrt{E} + C + N/E + E_{2}/E^{2} + D\ln E$",
}


# ---------------------------------------------------------------------------
# Data access
# ---------------------------------------------------------------------------


def default_json_path():
    hist = os.environ.get("PATH_TO_HISTOGRAMS")
    if hist:
        return os.path.join(hist, "plots", "dashboard_data.json")
    return os.environ.get("DASHBOARD_JSON", "dashboard_data.json")


def get_series(data, method, process):
    """Return the energy-JER entry for one method x process, or None."""
    mv = data["methods"].get(method)
    if not mv:
        return None
    pv = mv["processes"].get(process)
    if not pv:
        return None
    entry = (pv.get("energy") or {}).get("_all")
    if not entry or len(entry.get("mid_points", [])) < 2:
        return None
    return entry


def get_fit(entry, model):
    """Return the stored fit dict for `model`, or None if it isn't available."""
    fit = (entry.get("fits") or {}).get(model)
    if not fit or not fit.get("fit_x"):
        return None
    return fit


def base_model_of(fit, model):
    return fit.get("base_model") or model.replace("_weighted", "")


def format_params(fit, model):
    """'S=0.07 N=0.0 C=0.01 D=3.3e-04' in the slides' print order."""
    base = base_model_of(fit, model)
    labels = PARAM_DISPLAY.get(base)
    popt = fit.get("popt") or []
    if not labels or len(labels) != len(popt):
        return " ".join(f"p{i}={v:.3g}" for i, v in enumerate(popt))
    by_label = dict(zip(labels, popt))
    parts = []
    for name in PRINT_ORDER:
        if name not in by_label:
            continue
        v = by_label[name]
        parts.append(f"{name}={v:.1e}" if name in SCIENTIFIC else f"{name}={round(v, 2)}")
    return " ".join(parts)


SHORT_NAMES = {
    "three_param": "3-param",
    "logE_4param": "4-param +log",
    "invE15_noise": "4-param +$E^{-1.5}$",
    "invE2_noise": "4-param +$E^{-2}$",
    "logE_and_invE2": "5-param +log+$E^{-2}$",
    "logE_invE2_noStoch": "4-param no-stoch",
    "two_param": "2-param",
}


def short_model_name(model):
    """'4-param +log (w)' - compact enough for a per-panel legend."""
    base = model.replace("_weighted", "")
    name = SHORT_NAMES.get(base, base)
    return f"{name} (w)" if model.endswith("_weighted") else name


def model_title(model, fit=None):
    base = base_model_of(fit or {}, model)
    formula = MODEL_FORMULA.get(base, base)
    suffix = " (error-weighted)" if model.endswith("_weighted") else ""
    return f"{formula}{suffix}"


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------


def draw_panel(ax, entry, fit, model, color, label, linestyle="--", errorbars=False,
               show_fit=True, markersize=4):
    """Points (+ fitted curve) for one method on one axis. Returns max y of points.

    With show_fit=False no curve is drawn and the legend entry is attached to the
    markers themselves, so the panel shows the measured resolution only.
    """
    x = np.asarray(entry["mid_points"], dtype=float)
    y = np.asarray(entry["sigma_over_E"], dtype=float)
    point_label = None if show_fit else label
    if errorbars and entry.get("sigma_over_E_err"):
        err = np.asarray(entry["sigma_over_E_err"], dtype=float)
        err = np.where(np.isfinite(err) & (err > 0), err, 0.0)
        ax.errorbar(x, y, yerr=err, fmt="x", markersize=markersize, color=color,
                    elinewidth=0.8, capsize=1.5, label=point_label)
    else:
        ax.plot(x, y, "x", markersize=markersize, color=color, label=point_label)
    if not show_fit:
        return float(np.nanmax(y)) if y.size else 0.0
    if fit is not None:
        ax.plot(fit["fit_x"], fit["fit_y"], linestyle=linestyle, color=color,
                label=f"{label} {format_params(fit, model)}")
    else:
        ax.plot([], [], linestyle=linestyle, color=color, label=f"{label} (fit n/a)")
    return float(np.nanmax(y)) if y.size else 0.0


def style_jer_axis(ax, ymax=None, headroom=1.25):
    ax.set_xlabel("$E_{true}$ [GeV]")
    ax.set_ylabel(r"$\sigma_E / E_{true}$")
    ax.grid()
    if ymax is not None and np.isfinite(ymax) and ymax > 0:
        # extra headroom keeps the (wide) parameter legend off the curves
        ax.set_ylim(0, headroom * ymax)


def annotate_matrix_plot_with_arrows(fig):
    """The two grey guide arrows on the 5x3 'full plots' slides."""
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax.set_axis_off()
    ax.annotate("", xy=(0.009, 0.60), xycoords="figure fraction",
                xytext=(0.009, 0.99), textcoords="figure fraction",
                arrowprops=dict(arrowstyle="->", lw=1.2), color="gray")
    ax.text(0.0048, 0.75, "More B-hadron content", transform=ax.transAxes,
            rotation=90, ha="center", va="center", fontsize=9.5)
    ax.annotate("", xy=(0.40, 0.99), xycoords="figure fraction",
                xytext=(0.01, 0.99), textcoords="figure fraction",
                arrowprops=dict(arrowstyle="->", lw=1.2), color="gray")
    ax.text(0.25, 0.995, "Higher number of final-state jets",
            transform=ax.transAxes, ha="center", va="center", fontsize=9.5)


def grid_shape(meta):
    rows = max(m["grid_row"] for m in meta.values()) + 1
    cols = max(m["grid_col"] for m in meta.values()) + 1
    return rows, cols


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def fig_grid_twin_scales(data, left_method, right_method, x_label_methods=None,
                         log_y=True):
    """Points-only grid with the two methods on independent y axes.

    Left axis = `left_method`, right axis = `right_method`, each with its own
    range shared across the whole grid. Every panel then uses the full height for
    both methods, at the cost of the vertical gap between the two sets of points
    no longer meaning anything - read each against its own (colour-matched) axis.
    """
    meta = data["process_meta"]
    rows, cols = grid_shape(meta)
    fig, ax = plt.subplots(rows, cols, figsize=(13, 11))
    used = {}
    lims = {left_method: [np.inf, 0.0], right_method: [np.inf, 0.0]}

    for process, pm in meta.items():
        r, c = pm["grid_row"], pm["grid_col"]
        a_left = ax[r, c]
        a_right = a_left.twinx()
        drew = False
        for method, a in ((left_method, a_left), (right_method, a_right)):
            entry = get_series(data, method, process)
            if entry is None:
                continue
            color = METHOD_COLORS.get(method, "black")
            x = np.asarray(entry["mid_points"], dtype=float)
            y = np.asarray(entry["sigma_over_E"], dtype=float)
            a.plot(x, y, "x", markersize=6, color=color,
                   label=METHOD_LABELS.get(method, method))
            good = y[np.isfinite(y) & (y > 0)]
            if good.size:
                lims[method][0] = min(lims[method][0], float(good.min()))
                lims[method][1] = max(lims[method][1], float(good.max()))
            drew = True
        if not drew:
            a_right.set_axis_off()
            continue
        used[(r, c)] = (a_left, a_right)
        a_left.set_title(pm["label"], fontsize=11)

    for (r, c), (a_left, a_right) in used.items():
        for method, a, side in ((left_method, a_left, "left"),
                                (right_method, a_right, "right")):
            color = METHOD_COLORS.get(method, "black")
            lo, hi = lims[method]
            if log_y:
                a.set_yscale("log")
            if np.isfinite(lo) and lo > 0 and hi > 0:
                a.set_ylim((0.8 * lo, 4.0 * hi) if log_y else (0, 1.2 * hi))
            label = r"$\sigma_E / E_{true}$" + f" ({METHOD_LABELS.get(method, method)})"
            a.set_ylabel(label, color=color, fontsize=9)
            a.tick_params(axis="y", colors=color, labelsize=8)
        a_left.set_xlabel("$E_{true}$ [GeV]")
        a_left.grid(alpha=0.4)
        handles = a_left.get_legend_handles_labels()[0] + \
            a_right.get_legend_handles_labels()[0]
        labels = a_left.get_legend_handles_labels()[1] + \
            a_right.get_legend_handles_labels()[1]
        a_left.legend(handles, labels, fontsize=9, loc="upper right", framealpha=0.95)

    for r in range(rows):
        for c in range(cols):
            if (r, c) not in used:
                ax[r, c].set_axis_off()
    fig.tight_layout()
    # No guide arrows: the per-panel jet-count badge carries that
    # information, and the arrows only crowded the canvas.
    return fig


def fig_grid(data, methods, model, errorbars=False, shared_ylim=True, show_fit=True,
             log_y=False):
    """The 5x3 process grid ('full plots' slides), one curve per method.

    With show_fit=False only the measured points are drawn. The shared y range is
    then computed from the methods passed in, so calling this once per detector
    family (PF vs. Calo) gives each family its own standardized scale instead of
    letting the much worse CaloJets resolution squash the PF panels.

    log_y=True instead puts everything on one logarithmic scale, which is the way
    to show PF and Calo together: the ~8x gap between them at low E is legible on
    a log axis without either being flattened.
    """
    meta = data["process_meta"]
    rows, cols = grid_shape(meta)
    fig, ax = plt.subplots(rows, cols, figsize=(13, 11))
    used = set()
    ymax_global = 0.0
    ymin_global = np.inf
    panel_ymax = {}
    panel_ymin = {}
    for process, pm in meta.items():
        r, c = pm["grid_row"], pm["grid_col"]
        a = ax[r, c]
        local_max = 0.0
        local_min = np.inf
        drew = False
        for method in methods:
            entry = get_series(data, method, process)
            if entry is None:
                continue
            fit = get_fit(entry, model) if show_fit else None
            local_max = max(local_max, draw_panel(
                a, entry, fit, model, METHOD_COLORS.get(method, "black"),
                METHOD_LABELS.get(method, method), errorbars=errorbars,
                show_fit=show_fit, markersize=4 if show_fit else 6))
            y = np.asarray(entry["sigma_over_E"], dtype=float)
            positive = y[np.isfinite(y) & (y > 0)]
            if positive.size:
                local_min = min(local_min, float(positive.min()))
            drew = True
        if not drew:
            continue
        used.add((r, c))
        panel_ymax[(r, c)] = local_max
        panel_ymin[(r, c)] = local_min
        ymax_global = max(ymax_global, local_max)
        ymin_global = min(ymin_global, local_min)
        a.set_title(pm["label"], fontsize=10 if show_fit else 11)
        # Jet multiplicity in the corner, so a single panel or row is readable
        # on its own once it is lifted out of the grid.
        # Upper left: the legend owns the upper right, and the lower left is not
        # free either - on the shared CaloJets scale the PF curves sit there.
        a.text(0.035, 0.96, f"{pm.get('n_jets', '?')} jets", transform=a.transAxes,
               ha="left", va="top", fontsize=14.5, color="darkred",
               bbox=dict(facecolor="white", edgecolor="red", linewidth=1.2,
                         boxstyle="round,pad=0.3"))
        # Without the fit-parameter text the legend is short, so it can be bigger.
        # Upper right is the corner the falling resolution leaves free, on both
        # linear and log axes (on log, the PF points sit along the bottom, so
        # lower left is not free even though it looks it).
        a.legend(fontsize=6.0 if show_fit else 10.0, loc="upper right",
                 framealpha=0.95)

    # 1.5 of headroom is there to keep the wide fit-parameter legend off the
    # curves; points-only panels don't need nearly as much.
    headroom = 1.5 if show_fit else 1.18
    for r in range(rows):
        for c in range(cols):
            if (r, c) not in used:
                ax[r, c].set_axis_off()
                continue
            if log_y:
                style_jer_axis(ax[r, c])  # labels/grid only, no linear ylim
                ax[r, c].set_yscale("log")
                lo = ymin_global if shared_ylim else panel_ymin[(r, c)]
                hi = ymax_global if shared_ylim else panel_ymax[(r, c)]
                if np.isfinite(lo) and lo > 0 and hi > 0:
                    # Extra top room so the legend clears the low-E Calo points.
                    ax[r, c].set_ylim(0.8 * lo, 4.0 * hi)
                continue
            style_jer_axis(ax[r, c],
                           ymax_global if shared_ylim else panel_ymax[(r, c)],
                           headroom=headroom)
    fig.tight_layout()
    # No guide arrows: the per-panel jet-count badge carries that
    # information, and the arrows only crowded the canvas.
    return fig


def fig_selected(data, methods, model, processes, errorbars=False):
    """The two stacked panels from the 'Jet Energy Resolution (Durham PF)' slide."""
    meta = data["process_meta"]
    fig, ax = plt.subplots(len(processes), 1, figsize=(5.2, 3.1 * len(processes)),
                           squeeze=False)
    ymax = 0.0
    for i, process in enumerate(processes):
        a = ax[i, 0]
        for method in methods:
            entry = get_series(data, method, process)
            if entry is None:
                continue
            ymax = max(ymax, draw_panel(
                a, entry, get_fit(entry, model), model,
                METHOD_COLORS.get(method, "black"),
                METHOD_LABELS.get(method, method), errorbars=errorbars))
        a.set_title(meta[process]["label"], fontsize=10)
        a.legend(fontsize=6.5)
    for i in range(len(processes)):
        style_jer_axis(ax[i, 0], ymax)
    fig.tight_layout()
    return fig


def fig_overlay(data, method, model, processes, group_title, errorbars=False):
    """JER (left) + energy response (right) for a set of processes, one method.

    Reproduces the 'Final state containing b-jets' slide: processes overlaid in
    their own colours, with the response ratio next to it.
    """
    meta = data["process_meta"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ymax = 0.0
    for process in processes:
        entry = get_series(data, method, process)
        if entry is None:
            continue
        pm = meta[process]
        color = pm["color"]
        ymax = max(ymax, draw_panel(ax[0], entry, get_fit(entry, model), model,
                                    color, pm["label"], errorbars=errorbars))
        response = entry.get("response")
        if response:
            ax[1].plot(entry["mid_points"], response, "--x", color=color,
                       markersize=4, label=pm["label"])
    style_jer_axis(ax[0], ymax)
    ax[0].set_title(group_title, fontsize=11)
    ax[0].legend(fontsize=7, title=f"{METHOD_LABELS.get(method, method)}",
                 title_fontsize=7.5)
    ax[1].set_xlabel("$E_{true}$ [GeV]")
    ax[1].set_ylabel("Response in $E_{reco}$ / $E_{true}$")
    ax[1].grid()
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    return fig


def load_mass_hist(data, plots_dir, method, process, definition="reco", full=True):
    """(edges, y) for one mH histogram; prefers the full-resolution sidecar file.

    The dashboard stores a 200-bin downsample inline plus a `full_url` pointing
    at the 1000-bin version - the peak is only ~1.3 GeV wide, so the downsample
    visibly under-resolves it. y is already normalised to unit area.
    """
    pv = (data["methods"].get(method) or {}).get("processes", {}).get(process)
    entry = ((pv or {}).get("mass") or {}).get("definitions", {}).get(definition)
    if not entry:
        return None
    if full and entry.get("full_url"):
        path = os.path.join(plots_dir, entry["full_url"])
        if os.path.exists(path):
            with open(path) as fh:
                full_entry = json.load(fh)
            return (np.asarray(full_entry["edges"], dtype=float),
                    np.asarray(full_entry["y"], dtype=float))
    return (np.asarray(entry["edges"], dtype=float),
            np.asarray(entry["y"], dtype=float))


def fig_mass_overlay(data, plots_dir, method, processes, group_title,
                     legend_title=None, xlim=(80.0, 150.0), definition="reco",
                     full=True):
    """One single-panel mH figure for a set of processes.

    Replaces the slides' two-panel (full range + zoom) layout: the full-range
    panel is almost entirely empty, so only the zoom is kept.
    """
    meta = data["process_meta"]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ymax = 0.0
    for process in processes:
        hist = load_mass_hist(data, plots_dir, method, process, definition, full)
        if hist is None:
            continue
        edges, y = hist
        pm = meta[process]
        ax.stairs(y, edges, color=pm["color"], lw=1.4, label=pm["label"])
        centers = 0.5 * (edges[:-1] + edges[1:])
        inside = (centers >= xlim[0]) & (centers <= xlim[1])
        if inside.any():
            ymax = max(ymax, float(np.nanmax(y[inside])))
    ax.set_xlim(*xlim)
    if ymax > 0:
        # Autoscale to the tallest peak *inside* the x-window (plain autoscale
        # would also see the out-of-window part of the histogram).
        ax.set_ylim(0, 1.06 * ymax)
    ax.set_xlabel("$m_H$ [GeV]")
    ax.set_ylabel("Normalized events")
    ax.grid(alpha=0.5)
    ax.set_title(group_title, fontsize=12)
    # Peak sits right of centre with the tail to the left, so upper-left is the
    # one corner that is reliably empty for every process set here.
    ax.legend(fontsize=9.5, loc="upper left", title=legend_title,
              title_fontsize=9.5, framealpha=0.92, borderpad=0.5,
              labelspacing=0.35, handlelength=1.6)
    fig.tight_layout()
    return fig


def fig_model_comparison(data, method, model, compare_model, errorbars=False):
    """Same 5x3 grid, but overlaying two fit models for a single method.

    The backup-slide style check that the extra log term actually helps: the
    legend carries each model's chi2/ndf where the dashboard stored one.
    """
    meta = data["process_meta"]
    rows, cols = grid_shape(meta)
    fig, ax = plt.subplots(rows, cols, figsize=(13, 11))
    styles = {model: ("crimson", "-"), compare_model: ("dimgray", "--")}
    short = {model: short_model_name(model), compare_model: short_model_name(compare_model)}
    used, ymax_global, panel_ymax = set(), 0.0, {}
    for process, pm in meta.items():
        r, c = pm["grid_row"], pm["grid_col"]
        a = ax[r, c]
        entry = get_series(data, method, process)
        if entry is None:
            continue
        x = np.asarray(entry["mid_points"], dtype=float)
        y = np.asarray(entry["sigma_over_E"], dtype=float)
        if errorbars and entry.get("sigma_over_E_err"):
            err = np.asarray(entry["sigma_over_E_err"], dtype=float)
            err = np.where(np.isfinite(err) & (err > 0), err, 0.0)
            a.errorbar(x, y, yerr=err, fmt="x", markersize=4, color="black",
                       elinewidth=0.8, capsize=1.5)
        else:
            a.plot(x, y, "x", markersize=4, color="black")
        for m in (model, compare_model):
            fit = get_fit(entry, m)
            if fit is None:
                continue
            color, ls = styles[m]
            chi2 = fit.get("chi2_ndf")
            tag = f", $\\chi^2$/ndf={chi2:.3g}" if chi2 is not None else ""
            a.plot(fit["fit_x"], fit["fit_y"], linestyle=ls, color=color,
                   label=f"{short[m]}: {format_params(fit, m)}{tag}")
        used.add((r, c))
        panel_ymax[(r, c)] = float(np.nanmax(y)) if y.size else 0.0
        ymax_global = max(ymax_global, panel_ymax[(r, c)])
        a.set_title(pm["label"], fontsize=10)
        a.legend(fontsize=6.0, loc="upper right",
                 title=METHOD_LABELS.get(method, method), title_fontsize=6.5)
    for r in range(rows):
        for c in range(cols):
            if (r, c) not in used:
                ax[r, c].set_axis_off()
                continue
            style_jer_axis(ax[r, c], ymax_global, headroom=1.6)
    fig.tight_layout()
    # No guide arrows: the per-panel jet-count badge carries that
    # information, and the arrows only crowded the canvas.
    return fig


def fig_residuals(data, methods, model, compare_model):
    """Fit - data residuals vs E, all processes overlaid, one panel per method.

    Quantifies the slides' 'still working on the fits' caveat: how far each
    model sits from the measured points across the whole sample.
    """
    fig, ax = plt.subplots(1, len(methods), figsize=(5.2 * len(methods), 4.0),
                           squeeze=False)
    for i, method in enumerate(methods):
        a = ax[0, i]
        stats = {}
        for m, color in ((model, "crimson"), (compare_model, "dimgray")):
            resid_all = []
            for process in data["process_meta"]:
                entry = get_series(data, method, process)
                if entry is None:
                    continue
                fit = get_fit(entry, m)
                if fit is None:
                    continue
                x = np.asarray(entry["mid_points"], dtype=float)
                y = np.asarray(entry["sigma_over_E"], dtype=float)
                yhat = np.interp(x, fit["fit_x"], fit["fit_y"])
                resid = yhat - y
                a.plot(x, resid, "o", markersize=3, color=color, alpha=0.55)
                resid_all.append(resid)
            if resid_all:
                r = np.concatenate(resid_all)
                stats[m] = float(np.mean(np.abs(r)))
                a.plot([], [], "o", color=color,
                       label=f"{short_model_name(m)}  MAE={stats[m]:.2e}")
        a.axhline(0.0, color="black", lw=0.8)
        a.set_xlabel("$E_{true}$ [GeV]")
        a.set_ylabel(r"fit $-$ data  ($\sigma_E/E$)")
        a.set_title(METHOD_LABELS.get(method, method), fontsize=10)
        a.grid()
        a.legend(fontsize=7.5)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def processes_by_line_style(meta, style):
    return [p for p, m in meta.items() if m.get("line_style") == style]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", default=None,
                        help="dashboard_data.json (default: $PATH_TO_HISTOGRAMS/plots/dashboard_data.json)")
    parser.add_argument("--model", default=DEFAULT_MODEL,
                        help=f"fit model to draw (default: {DEFAULT_MODEL})")
    parser.add_argument("--compare-model", default="three_param",
                        help="second model for the comparison/residual figures")
    parser.add_argument("--outputDir", default=None,
                        help="output folder (default: <plots dir>/presentation_<model>)")
    parser.add_argument("--errorbars", action="store_true",
                        help="draw the stored per-point stat errors")
    parser.add_argument("--per-panel-ylim", action="store_true",
                        help="autoscale each grid panel instead of sharing one y range")
    parser.add_argument("--list-models", action="store_true",
                        help="print the fit models available in the JSON and exit")
    parser.add_argument("--points-only-dir", default=None,
                        help="folder for the fit-free JER grids "
                             "(default: <plots dir>/JER_points_only)")
    parser.add_argument("--mass-method", default="PF_Durham",
                        help="method used for the Higgs-mass figures (default: PF_Durham)")
    parser.add_argument("--mass-definition", default="reco",
                        choices=["reco", "gen", "gt", "gt_recomatched"],
                        help="which mH definition to draw (default: reco)")
    parser.add_argument("--mass-range", nargs=2, type=float, default=[80.0, 150.0],
                        metavar=("MIN", "MAX"), help="mH x-range (default: 80 150)")
    args = parser.parse_args()

    json_path = args.data or default_json_path()
    with open(json_path) as fh:
        data = json.load(fh)

    if args.list_models:
        found = {}
        for mv in data["methods"].values():
            for pv in mv["processes"].values():
                entry = (pv.get("energy") or {}).get("_all") or {}
                for name in (entry.get("fits") or {}):
                    found[name] = found.get(name, 0) + 1
        for name, n in sorted(found.items(), key=lambda kv: -kv[1]):
            print(f"{name:28s} {n} series")
        return

    outputDir = args.outputDir or os.path.join(os.path.dirname(os.path.abspath(json_path)),
                                               f"presentation_{args.model}")
    os.makedirs(outputDir, exist_ok=True)
    print(f"data:   {json_path}")
    print(f"model:  {args.model}")
    print(f"output: {outputDir}")

    meta = data["process_meta"]
    plots_dir = os.path.dirname(os.path.abspath(json_path))  # full_hist/ lives here
    pf_pair = ["PF_Durham", "PF_Durham_IdealMatching"]
    pf_calo = ["PF_Durham", "CaloJets_Durham", "PF_Durham_IdealMatching"]
    shared = not args.per_panel_ylim

    figures = [
        ("JER_grid_PFJets.pdf",
         lambda: fig_grid(data, ["PF_Durham"], args.model, args.errorbars, shared)),
        ("JER_grid_PFJets_CaloJets.pdf",
         lambda: fig_grid(data, ["PF_Durham", "CaloJets_Durham"], args.model,
                          args.errorbars, shared)),
        ("JER_grid_PFJets_vs_IdealMatching.pdf",
         lambda: fig_grid(data, pf_pair, args.model, args.errorbars, shared)),
        ("JER_grid_PFJets_CaloJets_IdealMatching.pdf",
         lambda: fig_grid(data, pf_calo, args.model, args.errorbars, shared)),
        ("JER_selected_2jet_processes.pdf",
         lambda: fig_selected(data, pf_pair, args.model, SELECTED_PROCESSES, args.errorbars)),
        ("JER_overlay_bjet_processes.pdf",
         lambda: fig_overlay(data, "PF_Durham", args.model,
                             processes_by_line_style(meta, "-"),
                             "Final state containing b-jets", args.errorbars)),
        ("JER_overlay_lightflavour_processes.pdf",
         lambda: fig_overlay(data, "PF_Durham", args.model,
                             processes_by_line_style(meta, ":"),
                             "Light-flavour final state (q $\\in$ {u, d, s})", args.errorbars)),
        ("JER_grid_fit_comparison_PFJets.pdf",
         lambda: fig_model_comparison(data, "PF_Durham", args.model,
                                      args.compare_model, args.errorbars)),
        ("JER_grid_fit_comparison_CaloJets.pdf",
         lambda: fig_model_comparison(data, "CaloJets_Durham", args.model,
                                      args.compare_model, args.errorbars)),
        ("JER_fit_residuals.pdf",
         lambda: fig_residuals(data, pf_calo, args.model, args.compare_model)),
        ("mH_bjet_processes.pdf",
         lambda: fig_mass_overlay(data, plots_dir, args.mass_method,
                                  processes_by_line_style(meta, "-"),
                                  "Final state containing b-jets",
                                  legend_title=METHOD_LABELS.get(args.mass_method,
                                                                 args.mass_method),
                                  xlim=tuple(args.mass_range),
                                  definition=args.mass_definition)),
        ("mH_lightflavour_processes.pdf",
         lambda: fig_mass_overlay(data, plots_dir, args.mass_method,
                                  processes_by_line_style(meta, ":"),
                                  "H $\\rightarrow$ light-flavour jets",
                                  legend_title=f"{METHOD_LABELS.get(args.mass_method, args.mass_method)}"
                                               "\nq $\\in$ {u, d, s}",
                                  xlim=tuple(args.mass_range),
                                  definition=args.mass_definition)),
    ]

    written = []
    for name, builder in figures:
        try:
            fig = builder()
        except Exception as exc:  # keep going so one bad figure can't kill the set
            print(f"  ⚠️ {name}: {type(exc).__name__}: {exc}")
            continue
        path = os.path.join(outputDir, name)
        fig.savefig(path)
        plt.close(fig)
        written.append((name, FIGURE_DESCRIPTIONS.get(name, "")))
        print(f"  wrote {path}")

    # Points-only grids, in their own folder: no fit is involved, so they don't
    # belong under presentation_<model>. Each detector family gets its own shared
    # y range (see fig_grid) - one scale for the PF grids, one for CaloJets,
    # standardized across all panels within each.
    points_dir = args.points_only_dir or os.path.join(
        os.path.dirname(os.path.abspath(json_path)), "JER_points_only")
    os.makedirs(points_dir, exist_ok=True)
    # errorbars=False throughout: these figures are meant to be the measured
    # points and nothing else (the stored stat errors are smaller than the
    # markers nearly everywhere anyway). --errorbars still applies to the fit
    # figures above.
    points_figures = [
        ("JER_grid_points_PFJets.pdf",
         lambda: fig_grid(data, ["PF_Durham"], args.model, errorbars=False,
                          shared_ylim=True, show_fit=False),
         "points only, PFJets"),
        ("JER_grid_points_PFJets_CaloJets.pdf",
         lambda: fig_grid(data, ["PF_Durham", "CaloJets_Durham"], args.model,
                          errorbars=False, shared_ylim=True, show_fit=False),
         "points only, PFJets + CaloJets, one shared linear y range"),
        ("JER_grid_points_PFJets_vs_IdealMatching.pdf",
         lambda: fig_grid(data, pf_pair, args.model, errorbars=False,
                          shared_ylim=True, show_fit=False),
         "points only, PFJets + PFJets/Ideal Matching, one shared y range"),
        ("JER_grid_points_CaloJets.pdf",
         lambda: fig_grid(data, ["CaloJets_Durham"], args.model, errorbars=False,
                          shared_ylim=True, show_fit=False),
         "points only, CaloJets, its own shared y range"),
        # Log y is only used for the PF+Calo figure: it exists precisely because
        # those two need a common scale, and a linear one can't hold both.
        ("JER_grid_points_PFJets_CaloJets_IdealMatching.pdf",
         lambda: fig_grid(data, pf_calo, args.model, errorbars=False,
                          shared_ylim=True, show_fit=False),
         "points only, PFJets + CaloJets + Ideal Matching, one shared linear y "
         "range (CaloJets dominate it - see the log version)"),
        ("JER_grid_points_PFJets_CaloJets_log.pdf",
         lambda: fig_grid(data, ["PF_Durham", "CaloJets_Durham"], args.model,
                          errorbars=False, shared_ylim=True, show_fit=False,
                          log_y=True),
         "points only, PFJets vs. CaloJets on one shared log y scale"),
        ("JER_grid_points_PFJets_CaloJets_log_separateLegend.pdf",
         lambda: fig_grid_twin_scales(data, "PF_Durham", "CaloJets_Durham"),
         "points only, PFJets on the left log axis and CaloJets on their own "
         "right log axis (each range shared across the grid)"),
    ]
    for scan_methods, scan_name, scan_desc in (
        (AK_METHODS, "JER_grid_points_Durham_vs_AntiKt.pdf",
         "points only, Durham vs. the generalized e+e- anti-kt radius scan"),
        (AK_ER_METHODS, "JER_grid_points_Durham_vs_AntiKt_Erecovery.pdf",
         "points only, Durham vs. the anti-kt radius scan with energy recovery"),
    ):
        if not any(m in data["methods"] for m in scan_methods):
            print(f"  (no anti-kt methods in this tree - skipping {scan_name})")
            continue
        points_figures.append((
            scan_name,
            lambda ms=scan_methods: fig_grid(data, ["PF_Durham"] + ms, args.model,
                                             errorbars=False, shared_ylim=True,
                                             show_fit=False),
            scan_desc))

    points_written = []
    for name, builder, desc in points_figures:
        try:
            fig = builder()
        except Exception as exc:
            print(f"  ⚠️ {name}: {type(exc).__name__}: {exc}")
            continue
        path = os.path.join(points_dir, name)
        fig.savefig(path)
        plt.close(fig)
        points_written.append((name, desc))
        print(f"  wrote {path}")

    with open(os.path.join(points_dir, "README.md"), "w") as fh:
        fh.write("# JER grids, measured points only (no fits)\n\n")
        fh.write("Same 5x3 process grid as the presentation figures, showing the\n")
        fh.write("per-energy-bin `std68` resolution points with no fitted curve, so\n")
        fh.write("nothing here depends on a choice of fit model.\n\n")
        fh.write("The y range is standardized across all panels *within* a figure but\n")
        fh.write("differs between them: CaloJets are far enough off the PF resolution\n")
        fh.write("that a single common scale would flatten the PF panels.\n\n")
        fh.write("Markers only - no fitted curve and no error bars (the stored\n")
        fh.write("per-point stat errors are smaller than the markers nearly everywhere).\n\n")
        fh.write("Regenerate with:\n\n    python3 src/plotting/presentation_jer_plots.py\n\n")
        fh.write("| file | figure |\n|---|---|\n")
        for name, desc in points_written:
            fh.write(f"| `{name}` | {desc} |\n")
    print(f"  wrote {os.path.join(points_dir, 'README.md')}")

    readme = os.path.join(outputDir, "README.md")
    with open(readme, "w") as fh:
        fh.write(f"# Presentation JER figures - `{args.model}` fit\n\n")
        fh.write(f"Energy-dependence model: {model_title(args.model)}\n\n")
        fh.write(f"Comparison model in the comparison/residual figures: "
                 f"`{args.compare_model}`.\n\n")
        fh.write(f"Regenerate with:\n\n    python3 src/plotting/presentation_jer_plots.py"
                 f" --model {args.model}\n\n")
        fh.write("| file | figure |\n|---|---|\n")
        for name, desc in written:
            fh.write(f"| `{name}` | {desc} |\n")
    print(f"  wrote {readme}")


if __name__ == "__main__":
    main()
