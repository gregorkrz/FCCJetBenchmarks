"""Decompose the reconstructed Higgs-mass peak into its contributions.

The production mH figures show *how wide* the peak is, not *why*. This script
overlays a nested ladder of mH definitions, one panel per process, in the same
3x5 grid as the JER "full plots" slides:

    Visible final-state     the visible final-state gen particles from the H. What it
    particles               lacks is the invisible loss (neutrinos from semileptonic
                            heavy-flavour decays); the H natural width is zero here.
    Jet def.                gen jets with gen momenta: the cost of partitioning the
                            event into jets with a perfect detector. Mixes QCD
                            radiation crossing between the H and Z systems
                            (irreducible) with the clustering convention (Durham,
                            exclusive N, the dR matching radius).
    Detector                each gen particle replaced by its RecoMCLink partner:
                            efficiency, |eta| < 2.56 acceptance, smearing. No jets, so
                            no assignment error can enter.
    Detector + Jet def.     both of the above - gen jets whose constituents are
                            replaced by their reco partners (PF_Durham_IdealMatching).
    ... + Reco-jet          the reco particles themselves are clustered into jets and
    clustering              the H jets found by dR matching (PF_Durham) - the standard
                            reco mH. Undefined when the matching fails, so unlike the
                            others it sits on a subset of the events.

All four come out of the histmaker (1000 bins, 0-250 GeV), so this script only
reads and draws them - see VARIANTS below for the provenance of each curve. The
histograms are read straight out of the histmaker ROOT files with uproot (no
ROOT/container needed, and unlike the dashboard JSON this keeps the raw entry
counts, which are needed for the kept-fraction column).

Usage:
    source env.sh
    python src/plotting/mh_decomposition_plots.py
    python src/plotting/mh_decomposition_plots.py --rebin 4
    python src/plotting/mh_decomposition_plots.py --inputDir <subset> --reco-hist h_mH_reco_fixed
"""
import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# The `src.` imports below need the repo root on the path; running the script by
# file name only puts src/plotting/ there.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.process_config import (
    HUMAN_READABLE_PROCESS_NAMES,
    NUMBER_OF_JETS,
    PROCESS_COLORS,
    PROCESS_TO_ROW_COL,
)
from src.plotting.resolution_methods import fit_peak

# ---------------------------------------------------------------------------
# The ladder. `method` is "reco" or "ideal", resolved to a directory by the CLI
# args; `hist` of None means "the jet-based reco mH histogram" (--reco-hist), so
# a re-run with the fixed H-jet mapping can be plotted by name.
# ---------------------------------------------------------------------------

VARIANTS = [
    dict(
        key="physics",
        label="Visible final-state particles",
        method="reco",
        hist="h_mH_stable_gt_particles",
        color="black",
        linestyle="-",
        lw=0.8,
        # A 0.25 GeV-wide delta at 125 GeV: letting it drive the linear y-axis
        # would flatten every other curve, so it is allowed to clip instead.
        drives_ylim=False,
        note="visible final-state gen particles from H (neutrinos removed)",
    ),
    dict(
        # Jet definition with no detector at all, so it can be read against the
        # "Detector" curve to see how the two effects combine.
        key="jet_def_only",
        label="Jet def.",
        method="reco",
        hist="h_mH_gen",
        color="tab:purple",
        linestyle="-",
        lw=1.6,
        # For 2-jet processes the two gen jets recover almost the whole visible
        # mass, so this is near delta-like and would flatten those panels; for
        # 4/6 jets it is broad and never the tallest curve anyway.
        drives_ylim=False,
        note="gen jets with gen momenta - no detector at all",
    ),
    dict(
        key="detector",
        label="Detector",
        method="reco",
        hist="h_mH_reco_particles_matched",
        color="tab:green",
        linestyle="-",
        lw=1.4,
        note="those gen particles -> their RecoMCLink reco partners (no jets)",
    ),
    dict(
        key="det_phys",
        label="Detector + Jet def.",
        method="ideal",
        hist=None,
        color="tab:orange",
        linestyle="-",
        # Drawn thick and before the next one, which nearly coincides with it, so
        # the thinner line on top leaves an orange halo where they agree.
        lw=2.6,
        note="gen jets, constituents replaced by reco partners (ideal matching)",
    ),
    dict(
        key="det_phys_pf",
        label="Detector + Jet def. + Reco-jet clustering",
        method="reco",
        hist=None,
        color="tab:blue",
        linestyle="-",
        lw=1.2,
        note="jets clustered from the reco particles, H jets by dR matching (PF jets)",
    ),
]

LADDERS = {
    # figure stem -> variant keys drawn in it
    "mH_decomposition_3": ["detector", "det_phys", "det_phys_pf"],
    "mH_decomposition_3plus1": ["physics", "detector", "det_phys", "det_phys_pf"],
    # Everything, including the detector-free gen-jet curve: lets the
    # jet-definition and detector contributions be compared side by side.
    # No invisible-loss curve here: it is a 1-bin spike, negligible next to the
    # others, and it only forced the y axis to be clipped.
    "mH_decomposition_all": ["detector", "jet_def_only", "det_phys",
                             "det_phys_pf"],
    # Same figure with the gen-level visible-particle reference added.
    "mH_decomposition_all_visible": ["physics", "detector", "jet_def_only",
                                     "det_phys", "det_phys_pf"],
}

# A second, deliberately blunter label set for the same four curves - the version
# to put on a slide or poster, where "Physics" stands for the gen-level jet
# definition and "PF" for clustering the reco particles.
# "Detector + Physics + PF" is left out on purpose: it sits practically on top
# of "Detector + Physics" (the widths differ by <5% everywhere, see
# mh_decomposition_summary.md), so on a slide it only thickens the orange line.
# The mH_decomposition_all* figures still carry it.
SIMPLE_LABELS = {
    "jet_def_only": "Physics",
    "detector": "Detector",
    "det_phys": "Detector + Physics",
}
# Uniform thin lines: the thick orange in VARIANTS only exists so the PF curve
# drawn on top of it stays visible, and the PF curve isn't in this set.
VARIANTS_SIMPLE = [dict(v, label=SIMPLE_LABELS[v["key"]], lw=1.3)
                   for v in VARIANTS if v["key"] in SIMPLE_LABELS]
SIMPLE_LADDERS = {"mH_decomposition_simple": list(SIMPLE_LABELS)}

# A third label scheme, on the same three PFlow-based histograms: every curve is
# built from PFlow objects and they differ only in how those objects are grouped
# into the Higgs candidate, so no gen-momentum curve appears here.
#   Detector           perfect grouping from MC ancestry - width is detector only
#   Physics            grouping from the truth (gen) jets, constituents replaced
#                      by their PFlow partners
#   Detector + Physics grouping from the reco-clustered jets, i.e. everything reco
PFLOW_LABELS = {
    "detector": "Detector",
    "det_phys": "Physics",
    "det_phys_pf": "Detector + Physics",
}
VARIANTS_PFLOW = [dict(v, label=PFLOW_LABELS[v["key"]], lw=1.3)
                  for v in VARIANTS if v["key"] in PFLOW_LABELS]
PFLOW_LADDERS = {"mH_decomposition_pflow": list(PFLOW_LABELS)}

# Same three-way split as PFLOW_LABELS, but with the Physics curve taken from
# pure generator level (gen jets, gen momenta) instead of gen jets whose
# constituents were swapped for their PFlow partners. That makes the three
# curves one effect each rather than two of them sharing the detector:
#   Physics            jet definition with a perfect detector (no reco at all)
#   Detector           detector with perfect grouping (no jets at all)
#   Detector + Physics both together - fully reconstructed jets, the real mH
GENPHYS_LABELS = {
    "jet_def_only": "Physics",
    "detector": "Detector",
    "det_phys_pf": "Detector + Physics",
}
VARIANTS_GENPHYS = [dict(v, label=GENPHYS_LABELS[v["key"]], lw=1.3)
                    for v in VARIANTS if v["key"] in GENPHYS_LABELS]
GENPHYS_LADDERS = {"mH_decomposition_genphys": list(GENPHYS_LABELS)}

# Transposed layout: jet multiplicity along the row (2 -> 4 -> 6 jets), flavour
# down the rows, so a whole row can be lifted onto a poster as the "2, 4 and 6
# jets" comparison for one flavour. Only rows with a complete triple are useful
# for that, which is why these are listed explicitly rather than derived from
# PROCESS_TO_ROW_COL (whose b-rich column has no 6-jet entry).
JETS_IN_COLUMNS = [
    ["p8_ee_ZH_vvqq_ecm240", "p8_ee_ZH_qqqq_ecm240", "p8_ee_ZH_6jet_LF_ecm240"],
    ["p8_ee_ZH_vvgg_ecm240", "p8_ee_ZH_qqgg_ecm240", "p8_ee_ZH_6jet_ecm240"],
    ["p8_ee_ZH_vvbb_ecm240", "p8_ee_ZH_bbbb_ecm240", "p8_ee_ZH_6jet_HF_ecm240"],
]


def transposed_layout(data, jets_rows=None):
    """process -> (row, col) with 2/4/6 jets in the columns."""
    layout = {}
    for r, row in enumerate(jets_rows or JETS_IN_COLUMNS):
        for c, process in enumerate(row):
            if data.get(process):
                layout[process] = (r, c)
    return layout


# Per-variant single-canvas figures overlaying the three jet multiplicities. The
# two sets are the columns of the grid: the light-flavour one is its first column
# (grid_col 0), the b-rich one the rightmost occupied cell of each row - i.e. the
# two ends of the "more B-hadron content" arrow.
JET_MULTIPLICITY_SETS = {
    "lightflavour": dict(
        label="Light-flavour final states",
        processes=[
            ("2 jets", "p8_ee_ZH_vvqq_ecm240"),
            ("4 jets", "p8_ee_ZH_qqqq_ecm240"),
            ("6 jets", "p8_ee_ZH_6jet_LF_ecm240"),
        ],
    ),
    "bjetrich": dict(
        label="b-jet-rich final states",
        processes=[
            ("2 jets", "p8_ee_ZH_vvbb_ecm240"),
            ("4 jets", "p8_ee_ZH_bbbb_ecm240"),
            ("6 jets", "p8_ee_ZH_6jet_HF_ecm240"),
        ],
    ),
}

# The same three PFlow definitions, but restricted to events where *both* jet
# definitions found all the Higgs jets, so the three curves sit on one event
# sample and are directly comparable. All three histograms come out of a single
# PF_Durham run (see histmaker.py), which is the only way to intersect the
# selections - separate per-method histograms have already integrated over
# events and cannot be intersected after the fact.
# The genphys split (one effect per curve) on the common event sample.
GENPHYS_COMMON_VARIANTS = [
    dict(key="common_gen", label="Physics", method="reco",
         hist="h_mH_common_gen", color="tab:purple", linestyle="-", lw=1.3,
         drives_ylim=False,
         note="truth (gen) jets with gen momenta - no detector at all"),
    dict(key="common_detector", label="Detector", method="reco",
         hist="h_mH_common_detector", color="tab:green", linestyle="-", lw=1.3,
         note="perfect clustering (MC ancestry), PFlow objects"),
    dict(key="common_reco", label="Detector + Physics", method="reco",
         hist="h_mH_common_reco", color="tab:blue", linestyle="-", lw=1.3,
         note="reco-clustered jets + gen-reco jet matching (everything reco)"),
]
GENPHYS_COMMON_LADDERS = {
    "mH_decomposition_genphys_common": [v["key"] for v in GENPHYS_COMMON_VARIANTS]
}

PFLOW_COMMON_VARIANTS = [
    dict(key="common_detector", label="Detector", method="reco",
         hist="h_mH_common_detector", color="tab:green", linestyle="-", lw=1.3,
         note="perfect clustering (MC ancestry), PFlow objects"),
    dict(key="common_ideal", label="Physics", method="reco",
         hist="h_mH_common_ideal", color="tab:orange", linestyle="-", lw=1.3,
         note="truth (gen) jets, constituents -> their PFlow partners"),
    dict(key="common_reco", label="Detector + Physics", method="reco",
         hist="h_mH_common_reco", color="tab:blue", linestyle="-", lw=1.3,
         note="reco-clustered jets + gen-reco jet matching (everything reco)"),
]
PFLOW_COMMON_LADDERS = {
    "mH_decomposition_pflow_common": [v["key"] for v in PFLOW_COMMON_VARIANTS]
}

# --compare-fixed: the same two jet-based curves with the old and the corrected
# parton -> reco-jet mapping, to size up the indexing bug on identical events.
FIX_COMPARISON_VARIANTS = [
    dict(key="ideal_old", label="Ideal matching, old map", method="ideal",
         hist="h_mH_reco", color="tab:orange", linestyle="-", lw=2.6,
         note="h_mH_reco: parton->genjet composed with the reco->gen map"),
    dict(key="ideal_fixed", label="Ideal matching, fixed map", method="ideal",
         hist="h_mH_reco_fixed", color="darkred", linestyle="--", lw=1.2,
         note="h_mH_reco_fixed: parton->genjet composed with a true gen->reco map"),
    dict(key="pf_old", label="PF jets, old map", method="reco",
         hist="h_mH_reco", color="tab:blue", linestyle="-", lw=2.6,
         note="h_mH_reco: parton->genjet composed with the reco->gen map"),
    dict(key="pf_fixed", label="PF jets, fixed map", method="reco",
         hist="h_mH_reco_fixed", color="navy", linestyle="--", lw=1.2,
         note="h_mH_reco_fixed: parton->genjet composed with a true gen->reco map"),
]


# ---------------------------------------------------------------------------
# Data access
# ---------------------------------------------------------------------------


def default_input_dir():
    hist = os.environ.get("PATH_TO_HISTOGRAMS")
    if not hist:
        raise SystemExit(
            "PATH_TO_HISTOGRAMS is not set and --inputDir was not given (source env.sh)"
        )
    return hist


def rebin_hist(y, edges, factor):
    """Merge `factor` adjacent bins (dropping any remainder at the top)."""
    factor = int(factor)
    if factor <= 1:
        return np.asarray(y, dtype=float), np.asarray(edges, dtype=float)
    n_new = len(y) // factor
    y_new = np.asarray(y[: n_new * factor], dtype=float).reshape(n_new, factor).sum(axis=1)
    edges_new = np.asarray(edges, dtype=float)[:: factor][: n_new + 1]
    return y_new, edges_new


def cdf_stats(edges, y):
    """Median / q16 / q84 / 68% half-width straight from the CDF.

    More honest than a Gaussian sigma for these asymmetric shapes (the low tail
    from invisible losses, the high tail from Z-jet contamination).
    `cumsum(y)[i]` is the CDF at the *upper* edge of bin i, so it is interpolated
    against `edges[1:]`; using the bin centers instead would shift every quantile
    down by half a bin.
    """
    total = float(np.sum(y))
    if total <= 0:
        return None
    cdf = np.cumsum(y) / total
    q16, med, q84 = (float(np.interp(f, cdf, edges[1:])) for f in (0.16, 0.5, 0.84))
    return dict(median=med, q16=q16, q84=q84, hw68=0.5 * (q84 - q16))


def load_variant(root_path, hist_name):
    """Read one mH histogram: full-resolution counts, edges and entry counts.

    The underflow bin holds the events where `invariant_mass` returned -1, i.e.
    where the expected number of Higgs jets was not found (functions.h:1003), so
    `kept_frac` is the fraction of the file's event sample for which the
    observable is defined at all.
    """
    if not os.path.isfile(root_path):
        return None, f"missing file {root_path}"
    with uproot.open(root_path) as f:
        if hist_name not in [k.split(";")[0] for k in f.keys(cycle=False)]:
            return None, f"missing histogram {hist_name} in {os.path.basename(root_path)}"
        h = f[hist_name]
        y = h.values().astype(float)
        edges = h.axis().edges().astype(float)
        total = float(h.values(flow=True).sum())
    in_range = float(y.sum())
    if in_range <= 0:
        return None, f"empty histogram {hist_name} in {os.path.basename(root_path)}"
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = float(edges[1] - edges[0])
    entry = dict(
        y=y,
        edges=edges,
        centers=centers,
        bin_width=bin_width,
        entries=in_range,
        total_entries=total,
        kept_frac=in_range / total if total > 0 else float("nan"),
    )
    entry.update(cdf_stats(edges, y) or {})
    # Gaussian mu/sigma on the unit-area density, same fit as mass_plots.py
    density = y / in_range / bin_width
    fit = fit_peak(centers, density, model="gaussian")
    entry["mu"], entry["sigma"] = (float(fit[2][0]), float(fit[2][1])) if fit else (None, None)
    return entry, None


def collect(input_dir, method_dirs, reco_hist, processes, rebin, variants=VARIANTS):
    """{process: {variant_key: entry}}, plus a list of warnings."""
    data, warnings = {}, []
    for process in processes:
        per_variant = {}
        for variant in variants:
            method = method_dirs[variant["method"]]
            hist_name = variant["hist"] or reco_hist
            path = os.path.join(input_dir, method, f"{process}.root")
            entry, err = load_variant(path, hist_name)
            if entry is None:
                warnings.append(f"{process} / {variant['label']}: {err}")
                continue
            entry["y_plot"], entry["edges_plot"] = rebin_hist(entry["y"], entry["edges"], rebin)
            bin_width = entry["edges_plot"][1] - entry["edges_plot"][0]
            entry["density_plot"] = entry["y_plot"] / entry["entries"] / bin_width
            entry["method"] = method
            entry["hist_name"] = hist_name
            per_variant[variant["key"]] = entry
        if per_variant:
            data[process] = per_variant
        else:
            warnings.append(f"{process}: no variant could be loaded, panel left empty")
    return data, warnings


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------


def legend_label(variant, entry):
    """Just the curve name - the numbers live in mh_decomposition_summary.md."""
    return variant["label"]


def definitions_text(method_dirs, reco_hist, keys, variants=VARIANTS):
    lines = [r"$\bf{Curves}$ (histmaker output):"]
    for variant in variants:
        if variant["key"] not in keys:
            continue
        lines.append(
            f"{variant['label']}\n    {method_dirs[variant['method']]} / "
            f"{variant['hist'] or reco_hist}\n    {variant['note']}"
        )
    lines.append(
        "The curves are cumulative: each adds the named effect\n"
        "to the one above it. Each is normalized to unit area.\n"
        "The two jet curves are only filled when all H jets\n"
        "were matched, so they sit on a subset of the events.\n"
        r"Peak positions and widths: $\tt{mh\_decomposition\_summary.md}$"
    )
    if "physics" in keys:
        lines.append(
            "The invisible-loss peak is ~1 bin (0.25 GeV) wide and is\n"
            "clipped in the linear panels; see the log-y version."
        )
    return "\n".join(lines)


def annotate_grid_arrows(fig, transposed=False):
    """Guide arrows for the process grid, labelled to match the actual layout.

    Default layout: rows = jet multiplicity, columns = flavour. Transposed:
    columns = jet multiplicity, rows = flavour.
    (`presentation_jer_plots.annotate_matrix_plot_with_arrows()` hardcodes the
    transposed labels, which is why the JER "full plots" grids - built on
    `PROCESS_TO_ROW_COL`, i.e. the default layout - carry swapped ones.)
    """
    down_label = ("More B-hadron content" if transposed
                  else "Higher number of final-state jets")
    across_label = ("Higher number of final-state jets" if transposed
                    else "More B-hadron content")
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax.set_axis_off()
    ax.annotate("", xy=(0.008, 0.02), xycoords="figure fraction",
                xytext=(0.008, 0.955), textcoords="figure fraction",
                arrowprops=dict(arrowstyle="->", lw=1.2, color="gray"))
    ax.text(0.0038, 0.5, down_label,
            transform=ax.transAxes, rotation=90, ha="center", va="center",
            fontsize=10)
    ax.annotate("", xy=(0.55, 0.996), xycoords="figure fraction",
                xytext=(0.03, 0.996), textcoords="figure fraction",
                arrowprops=dict(arrowstyle="->", lw=1.2, color="gray"))
    # Label below the arrow line: above it there is no room left on the canvas.
    ax.text(0.29, 0.991, across_label, transform=ax.transAxes,
            ha="center", va="top", fontsize=10)


def packed_layout(data):
    """process -> (row, col), rows filled from the left with no gaps.

    Rows are still the jet multiplicity and the order within a row is still
    increasing B-hadron content (the `PROCESS_TO_ROW_COL` column), but the
    columns are packed so the only empty cells are trailing ones - a row can be
    cropped out of the figure for a slide without holes in the middle.
    """
    by_row = {}
    for process in data:
        if process not in PROCESS_TO_ROW_COL:
            continue
        r, c = PROCESS_TO_ROW_COL[process]
        by_row.setdefault(r, []).append((c, process))
    layout = {}
    for new_r, r in enumerate(sorted(by_row)):
        for new_c, (_, process) in enumerate(sorted(by_row[r])):
            layout[process] = (new_r, new_c)
    return layout


def fig_grid(data, keys, method_dirs, reco_hist, x_range, log_y, all_variants=VARIANTS,
             transposed=False, arrows=True, jets_rows=None, only_processes=None,
             jet_badge=False):
    """arrows=False drops the two grey guide arrows and the margin they need, so
    the panels run to the edge of the canvas - the version to drop on a poster.
    jets_rows overrides which rows of JETS_IN_COLUMNS are drawn (transposed only);
    only_processes restricts the grid to a subset of processes (either layout).
    jet_badge puts the jet multiplicity in the top-right corner of every panel,
    which is what makes a single row or column readable on its own.
    """
    loaded = {k for per_variant in data.values() for k in per_variant}
    keys = [k for k in keys if k in loaded]
    variants = [v for v in all_variants if v["key"] in keys]
    present = {p: v for p, v in data.items()
               if v and (only_processes is None or p in only_processes)}
    layout = (transposed_layout(present, jets_rows) if transposed
              else packed_layout(present))
    rows = max(rc[0] for rc in layout.values()) + 1
    cols = max(rc[1] for rc in layout.values()) + 1
    # Wider than tall per panel, like the mH comparison figures in joint_plots.
    fig, ax = plt.subplots(rows, cols, figsize=(4.0 * cols, 3.4 * rows),
                           squeeze=False)
    used = set()

    for process, rc in layout.items():
        per_variant = data.get(process)
        if not per_variant:
            continue
        r, c = rc
        a = ax[r, c]
        used.add(rc)
        # Autoscale from the bins inside the window only (plain autoscale would
        # also see the out-of-window part), and only from the curves flagged as
        # driving the y limit, so the delta-like invisible-loss spike doesn't
        # flatten everything else.
        ymax, ymax_all, ymin = 0.0, 0.0, np.inf
        for variant in variants:
            entry = per_variant.get(variant["key"])
            if entry is None:
                continue
            a.stairs(
                entry["density_plot"],
                entry["edges_plot"],
                color=variant["color"],
                linestyle=variant["linestyle"],
                lw=variant.get("lw", 1.3),
                label=legend_label(variant, entry),
            )
            centers = 0.5 * (entry["edges_plot"][:-1] + entry["edges_plot"][1:])
            inside = (centers >= x_range[0]) & (centers <= x_range[1])
            if not inside.any():
                continue
            vals = entry["density_plot"][inside]
            local_max = float(np.nanmax(vals))
            ymax_all = max(ymax_all, local_max)
            if variant.get("drives_ylim", True):
                ymax = max(ymax, local_max)
            positive = vals[vals > 0]
            if positive.size:
                ymin = min(ymin, float(np.nanmin(positive)))
        if ymax <= 0:
            ymax = ymax_all
        a.set_xlim(*x_range)
        if log_y and ymax_all > 0 and np.isfinite(ymin):
            a.set_yscale("log")
            a.set_ylim(max(ymin, ymax_all * 1e-4), ymax_all * 3.0)
        elif ymax > 0:
            a.set_ylim(0, 1.2 * ymax)
        a.set_title(HUMAN_READABLE_PROCESS_NAMES.get(process, process), fontsize=11)
        a.set_xlabel("$m_H$ [GeV]", fontsize=11)
        a.set_ylabel("Normalized events / GeV", fontsize=11)
        a.tick_params(axis="both", labelsize=9.5)
        a.grid(alpha=0.4)
        if jet_badge:
            a.text(0.97, 0.95, f"{NUMBER_OF_JETS.get(process, '?')} jets",
                   transform=a.transAxes, ha="right", va="top", fontsize=16,
                   color="darkred",
                   bbox=dict(facecolor="white", edgecolor="red", linewidth=1.2,
                             boxstyle="round,pad=0.35"))

    for r in range(rows):
        for c in range(cols):
            if (r, c) not in used:
                ax[r, c].set_axis_off()

    # One horizontal legend across the top for the whole grid, instead of a copy
    # in every panel. Narrow (transposed) figures wrap it onto a second line
    # rather than shrinking the text.
    handles = [Line2D([], [], color=v["color"], linestyle=v["linestyle"],
                      lw=max(v.get("lw", 1.3), 2.0), label=v["label"])
               for v in variants]
    if arrows:
        fig.tight_layout(rect=(0.022, 0, 1, 0.915))
        legend_y = 0.955
    else:
        # No arrow labels to clear: only the legend strip is reserved.
        legend_top = 1.0 - 0.55 / (3.4 * rows)
        fig.tight_layout(rect=(0, 0, 1, legend_top))
        legend_y = 1.0
    fig.legend(handles=handles, loc="upper center",
               ncol=len(handles) if cols >= 4 else 3,
               bbox_to_anchor=(0.5, legend_y), fontsize=11, frameon=True,
               framealpha=0.95, handlelength=2.2, columnspacing=2.0, borderpad=0.6)
    if arrows:
        annotate_grid_arrows(fig, transposed=transposed)
    return fig


def fig_decomposition_single(data, keys, process, x_range, log_y,
                             all_variants=VARIANTS):
    """One process, the whole ladder, on a single poster-sized canvas.

    Font sizes are set for a printed poster rather than for a 15-panel grid, so
    this is the figure to lift into a talk or poster for a given multiplicity.
    """
    per_variant = data.get(process)
    if not per_variant:
        return None
    variants = [v for v in all_variants if v["key"] in keys and v["key"] in per_variant]
    if not variants:
        return None
    fig, a = plt.subplots(figsize=(8.2, 6.0))
    ymax, ymax_all, ymin = 0.0, 0.0, np.inf
    for variant in variants:
        entry = per_variant[variant["key"]]
        a.stairs(entry["density_plot"], entry["edges_plot"], color=variant["color"],
                 linestyle=variant["linestyle"],
                 lw=max(variant.get("lw", 1.3) * 1.4, 1.8), label=variant["label"])
        centers = 0.5 * (entry["edges_plot"][:-1] + entry["edges_plot"][1:])
        inside = (centers >= x_range[0]) & (centers <= x_range[1])
        if not inside.any():
            continue
        vals = entry["density_plot"][inside]
        local_max = float(np.nanmax(vals))
        ymax_all = max(ymax_all, local_max)
        if variant.get("drives_ylim", True):
            ymax = max(ymax, local_max)
        positive = vals[vals > 0]
        if positive.size:
            ymin = min(ymin, float(np.nanmin(positive)))
    if ymax <= 0:
        ymax = ymax_all
    a.set_xlim(*x_range)
    if log_y and ymax_all > 0 and np.isfinite(ymin):
        a.set_yscale("log")
        a.set_ylim(max(ymin, ymax_all * 1e-4), ymax_all * 3.0)
    elif ymax > 0:
        a.set_ylim(0, 1.45 * ymax)
    a.set_xlabel("$m_H$ [GeV]", fontsize=17)
    a.set_ylabel("Normalized events / GeV", fontsize=17)
    a.tick_params(axis="both", labelsize=14)
    a.grid(alpha=0.45)
    a.set_title(
        f"{HUMAN_READABLE_PROCESS_NAMES.get(process, process)}"
        f"  -  {NUMBER_OF_JETS.get(process, '?')} jets", fontsize=19)
    a.legend(fontsize=14, loc="upper left", framealpha=0.94, labelspacing=0.5,
             handlelength=2.0, borderpad=0.7)
    fig.tight_layout()
    return fig


def write_per_process_pdf(data, output_dir, keys, views,
                          name="mH_decomposition_per_process.pdf"):
    """One poster-sized canvas per process, ordered 2-jet -> 4-jet -> 6-jet."""
    order = sorted(data, key=lambda p: (NUMBER_OF_JETS.get(p, 99),
                                        PROCESS_TO_ROW_COL.get(p, (9, 9))[1]))
    path = os.path.join(output_dir, name)
    pages = 0
    with PdfPages(path) as pdf:
        for process in order:
            for x_range, log_y in views:
                fig = fig_decomposition_single(data, keys, process,
                                               tuple(x_range), log_y)
                if fig is None:
                    continue
                scale = "log y" if log_y else "linear"
                fig.axes[0].set_title(
                    f"{HUMAN_READABLE_PROCESS_NAMES.get(process, process)}"
                    f"  -  {NUMBER_OF_JETS.get(process, '?')} jets"
                    f"  ({x_range[0]:g}-{x_range[1]:g} GeV, {scale})", fontsize=17)
                pdf.savefig(fig)
                plt.close(fig)
                pages += 1
    return path, pages


def fig_multiplicity(data, variant, process_set, x_range, log_y):
    """One canvas: a single mH definition, the 2/4/6-jet processes overlaid."""
    fig, a = plt.subplots(figsize=(6.6, 4.6))
    ymax, ymin = 0.0, np.inf
    drawn = []
    for label, process in process_set["processes"]:
        entry = (data.get(process) or {}).get(variant["key"])
        if entry is None:
            continue
        a.stairs(entry["density_plot"], entry["edges_plot"],
                 color=PROCESS_COLORS.get(process, "black"), lw=1.5, label=label)
        drawn.append((label, process))
        centers = 0.5 * (entry["edges_plot"][:-1] + entry["edges_plot"][1:])
        inside = (centers >= x_range[0]) & (centers <= x_range[1])
        if not inside.any():
            continue
        vals = entry["density_plot"][inside]
        ymax = max(ymax, float(np.nanmax(vals)))
        positive = vals[vals > 0]
        if positive.size:
            ymin = min(ymin, float(np.nanmin(positive)))
    if not drawn:
        plt.close(fig)
        return None
    a.set_xlim(*x_range)
    if log_y and ymax > 0 and np.isfinite(ymin):
        a.set_yscale("log")
        a.set_ylim(max(ymin, ymax * 1e-4), ymax * 3.0)
    elif ymax > 0:
        a.set_ylim(0, 1.25 * ymax)
    a.set_xlabel("$m_H$ [GeV]")
    a.set_ylabel("Normalized events / GeV")
    a.grid(alpha=0.5)
    a.set_title(f"{variant['label']}\n{process_set['label']}", fontsize=12)
    a.legend(fontsize=10, loc="upper left", framealpha=0.92, borderpad=0.5)
    # Which process each multiplicity actually is - the legend only says "N jets".
    caption = "   ".join(
        f"{label}: {HUMAN_READABLE_PROCESS_NAMES.get(process, process)}"
        for label, process in drawn
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.text(0.5, 0.015, caption, ha="center", fontsize=7.5)
    return fig


def write_multiplicity_pdf(data, output_dir, views, name="mH_jet_multiplicity.pdf"):
    """The per-variant multiplicity canvases as pages of a single PDF.

    `views` is a list of (x_range, log_y) - one page per view per
    (process set, mH definition).
    """
    path = os.path.join(output_dir, name)
    pages = 0
    with PdfPages(path) as pdf:
        for process_set in JET_MULTIPLICITY_SETS.values():
            for variant in VARIANTS:
                for x_range, log_y in views:
                    fig = fig_multiplicity(data, variant, process_set,
                                           tuple(x_range), log_y)
                    if fig is None:
                        continue
                    scale = "log y" if log_y else "linear"
                    fig.axes[0].set_title(
                        f"{variant['label']}\n{process_set['label']}"
                        f"  ({x_range[0]:g}-{x_range[1]:g} GeV, {scale})", fontsize=11)
                    pdf.savefig(fig)
                    plt.close(fig)
                    pages += 1
    return path, pages


# ---------------------------------------------------------------------------
# Tables / README
# ---------------------------------------------------------------------------

TABLE_COLUMNS = [
    "process", "variant", "method", "histogram", "entries", "kept_frac",
    "mu", "sigma", "median", "q16", "q84", "hw68",
]


def table_rows(data, variants=VARIANTS):
    rows = []
    for process, per_variant in data.items():
        for variant in variants:
            entry = per_variant.get(variant["key"])
            if entry is None:
                continue
            rows.append({
                "process": process,
                "variant": variant["label"],
                "method": entry["method"],
                "histogram": entry["hist_name"],
                "entries": f"{entry['entries']:.0f}",
                "kept_frac": f"{entry['kept_frac']:.4f}",
                "mu": "n/a" if entry["mu"] is None else f"{entry['mu']:.3f}",
                "sigma": "n/a" if entry["sigma"] is None else f"{entry['sigma']:.3f}",
                "median": f"{entry['median']:.3f}",
                "q16": f"{entry['q16']:.3f}",
                "q84": f"{entry['q84']:.3f}",
                "hw68": f"{entry['hw68']:.3f}",
            })
    return rows


def write_tables(rows, output_dir, stem="mh_decomposition_summary", title="mH decomposition summary"):
    csv_path = os.path.join(output_dir, f"{stem}.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TABLE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    md_path = os.path.join(output_dir, f"{stem}.md")
    with open(md_path, "w") as fh:
        fh.write(f"# {title}\n\n")
        fh.write("All masses in GeV. `hw68` = (q84 - q16)/2 from the CDF; `mu`/`sigma` from the\n")
        fh.write("Gaussian peak fit (argmax +/- 30 GeV window); `kept_frac` = fraction of the\n")
        fh.write("file's event sample for which the observable is defined.\n\n")
        fh.write("Widths below the 0.25 GeV bin width only measure the binning, not a physical\n")
        fh.write("width (a distribution contained in one bin gives hw68 = 0.085); they are shown\n")
        fh.write("as `<0.25` in the figure legends.\n\n")
        fh.write("| " + " | ".join(TABLE_COLUMNS) + " |\n")
        fh.write("| " + " | ".join("---" for _ in TABLE_COLUMNS) + " |\n")
        for row in rows:
            fh.write("| " + " | ".join(str(row[c]) for c in TABLE_COLUMNS) + " |\n")
    return csv_path, md_path


README_TEMPLATE = """# mH decomposition figures

Generated by `src/plotting/mh_decomposition_plots.py` from `{input_dir}`
(reco method: `{reco_method}`, ideal-matching method: `{ideal_method}`, jet-based mH
histogram: `{reco_hist}`, rebin factor: {rebin}).

One panel per process, 3x5 grid (rows = 2/4/6 final-state jets, columns = increasing
B-hadron content), same convention as the JER "full plots" slides. Every curve is
normalized to unit area, so the shapes are comparable even though the denominators
differ. The legends carry only the curve names; peak positions, widths and the
fraction of the event sample for which each definition is defined (`kept_frac`) are in
`mh_decomposition_summary.md` / `.csv`.

## Curves

| curve | method dir | histogram | meaning |
|---|---|---|---|
{variant_table}

The curves build up the effects one at a time:

- `Visible final-state particles` is the gen-level reference: the visible final-state
  particles descending from the Higgs. What it is missing is the invisible loss -
  neutrinos from semileptonic heavy-flavour decays - since the H natural width is zero
  in these samples. It is also the only curve without the |eta| < 2.56 cut.
- `Jet def.` clusters those gen particles into gen jets and uses gen momenta, so it is
  the cost of partitioning the event into jets with a perfect detector. That cost is
  itself a mixture: QCD radiation crossing between the H and Z systems and hadronization
  across colour strings (irreducible), plus the clustering convention - Durham,
  exclusive N, the dR matching radius (not irreducible; the anti-kt radius scan in
  `$PATH_TO_HISTOGRAMS` bounds that part).
- `Detector` uses no jets at all - each gen particle is replaced by its `RecoMCLink`
  partner - so no jet-assignment error can enter. It is reconstruction efficiency,
  acceptance and smearing only.
- `Detector + Jet def.` is both of the above: gen jets whose constituents are replaced
  by their reco partners (ideal matching).
- `Detector + Jet def. + Reco-jet clustering` additionally clusters the reco particles
  themselves and identifies the H jets by dR matching - the standard reco mH. Note this
  bundles particle-level misassignment with the jet matching, which these histograms
  cannot separate, and the matching also loses events (see caveat 3).
## Figures

{figure_list}

`mh_decomposition_summary.md` / `.csv` hold mu, sigma, median, q16, q84, 68% half-width
and the kept fraction for every process x curve.

## Caveats

{caveats}
"""

MAPPING_CAVEAT_OLD = """1. **H-jet selection index bug** (`src/histmaker_tools/truth_matching.py`, `merge_mappings`
   in `src/histmaker_functions/functions.h`): the parton->genjet map is composed with the
   **reco->gen** matching as if it were gen->reco. The selection is correct only where
   that permutation satisfies pi(g) = pi^-1(g), i.e. for jets sitting in cycles of length
   <= 2 (identity and pure swaps are fine). A cycle of length >= 3 picks a wrong but
   same-size jet subset, so the event survives with a wrong mass, typically in the high
   tail. **The 2-jet (Z->vv) row is provably unaffected**; only the 4-jet and 6-jet rows
   can be. This affects the two jet-based curves only. Removed by re-running the histmaker
   (which now also fills `h_mH_reco_fixed`) and passing `--reco-hist h_mH_reco_fixed`."""

MAPPING_CAVEAT_FIXED = """1. The jet-based curves use `{reco_hist}`, i.e. the corrected parton -> reco-jet mapping
   (a genuine gen->reco greedy matching instead of reusing the reco->gen one), so the
   H-jet-selection indexing bug does not apply here. Run with `--compare-fixed` to see
   old vs. fixed side by side on identical events."""

SAMPLE_CAVEAT = """2. **The two jet-based curves come from different files.** Their event samples differ by
   up to {worst:.1%} (worst process: `{worst_process}`), because each method dir applies its
   own fully-matched filter. The two jet-free curves are both taken from the
   `{reco_method}` file so that they at least share its sample. Re-run the two methods
   with `--no-filter-fully-matched` to remove this difference."""

SAMPLE_CAVEAT_OK = """2. The two jet-based curves come from different files, but their event samples agree to
   {worst:.2%}, so the comparison is effectively apples-to-apples. The two jet-free curves
   are both taken from the `{reco_method}` file."""


def sample_mismatch(data):
    """Worst relative disagreement in sample size between the two jet-based curves."""
    worst, worst_process = 0.0, "n/a"
    for process, per_variant in data.items():
        a = per_variant.get("det_phys")
        b = per_variant.get("det_phys_pf")
        if not a or not b or not b["total_entries"]:
            continue
        rel = abs(a["total_entries"] - b["total_entries"]) / b["total_entries"]
        if rel > worst:
            worst, worst_process = rel, process
    return worst, worst_process


def write_readme(output_dir, args, method_dirs, figures, data):
    variant_table = "\n".join(
        f"| {v['label']} | `{method_dirs[v['method']]}` | `{v['hist'] or args.reco_hist}` | {v['note']} |"
        for v in VARIANTS
    )
    figure_list = "\n".join(f"- `{name}` - {desc}" for name, desc in figures)
    worst, worst_process = sample_mismatch(data)
    caveats = "\n".join([
        (MAPPING_CAVEAT_FIXED if "fixed" in args.reco_hist else MAPPING_CAVEAT_OLD).format(
            reco_hist=args.reco_hist),
        (SAMPLE_CAVEAT_OK if worst < 0.01 else SAMPLE_CAVEAT).format(
            worst=worst, worst_process=worst_process, reco_method=args.reco_method),
        "3. The jet-based curves are only filled for events in which all Higgs jets were\n"
        "   parton-matched; `kept_frac` in the summary table is that fraction. Since every curve\n"
        "   is normalized to unit area, an apparent narrowing when going up the ladder can be\n"
        "   this selection rather than a genuine improvement.",
    ])
    with open(os.path.join(output_dir, "README.md"), "w") as fh:
        fh.write(README_TEMPLATE.format(
            input_dir=args.inputDir,
            reco_method=args.reco_method,
            ideal_method=args.ideal_method,
            reco_hist=args.reco_hist,
            rebin=args.rebin,
            variant_table=variant_table,
            figure_list=figure_list,
            caveats=caveats,
        ))


# ---------------------------------------------------------------------------
# Old vs. fixed H-jet mapping (--compare-fixed)
# ---------------------------------------------------------------------------


def moved_fraction(old, fixed):
    """Fraction of entries whose mH changed between the two mappings.

    Both histograms are filled on the same events, so sum|old - fixed| / 2 counts
    the events whose mass moved between bins. The two maps can also disagree on
    whether mH is defined at all (they hit `-1` for different events: the old map
    on unmatched *reco* jets, the fixed one on unmatched *gen* jets), which shows
    up as a difference in entry count rather than a difference in shape - without
    the extra term below such an event would only be counted as half a move.
    """
    if old is None or fixed is None or old["entries"] <= 0:
        return float("nan")
    moved = float(np.abs(old["y"] - fixed["y"]).sum())
    moved += abs(old["entries"] - fixed["entries"])
    return moved / 2.0 / old["entries"]


def fix_comparison(input_dir, method_dirs, args, processes, output_dir):
    """Extra figure + table: h_mH_reco vs h_mH_reco_fixed, same events."""
    print("compare-fixed: h_mH_reco vs h_mH_reco_fixed")
    data, warnings = collect(input_dir, method_dirs, args.reco_hist, processes,
                             args.rebin, variants=FIX_COMPARISON_VARIANTS)
    for w in warnings:
        print(f"  !! {w}")
    if not data:
        print("  !! no h_mH_reco_fixed histograms found, skipping the comparison")
        return []

    figures = []
    keys = [v["key"] for v in FIX_COMPARISON_VARIANTS]
    for suffix, x_range, log_y in (
        ("", args.zoom_range, False),
        ("_wide_log", args.wide_range, True),
    ):
        name = f"mH_mapping_fix_comparison{suffix}.pdf"
        fig = fig_grid(data, keys, method_dirs, args.reco_hist, tuple(x_range), log_y,
                       all_variants=FIX_COMPARISON_VARIANTS)
        fig.savefig(os.path.join(output_dir, name))
        plt.close(fig)
        scale = "log y" if log_y else "linear"
        figures.append((name, f"old vs fixed H-jet mapping, {x_range[0]:g}-{x_range[1]:g} GeV, {scale}"))
        print(f"  wrote {os.path.join(output_dir, name)}")

    rows = table_rows(data, variants=FIX_COMPARISON_VARIANTS)
    for path in write_tables(rows, output_dir, stem="mh_mapping_fix_summary",
                             title="Old vs. fixed parton -> reco-jet mapping"):
        print(f"  wrote {path}")

    print("  events whose mH changed with the fixed mapping (entries old -> fixed):")
    for process in processes:
        per_variant = data.get(process, {})
        parts = []
        for name, old_key, fixed_key in (("PF jets", "pf_old", "pf_fixed"),
                                         ("ideal", "ideal_old", "ideal_fixed")):
            old, fixed = per_variant.get(old_key), per_variant.get(fixed_key)
            frac = moved_fraction(old, fixed)
            counts = (f"{old['entries']:.0f}->{fixed['entries']:.0f}"
                      if old and fixed else "n/a")
            parts.append(f"{name} {frac:7.2%} ({counts})")
        print(f"    {process:26s} " + "   ".join(parts))
    return figures


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--inputDir", default=None,
                        help="histogram tree with one subdir per method (default: $PATH_TO_HISTOGRAMS)")
    parser.add_argument("--outputDir", default=None,
                        help="default: <inputDir>/plots/mh_decomposition")
    parser.add_argument("--reco-method", default="PF_Durham",
                        help="method dir for the PF-jet curve and the two jet-free curves")
    parser.add_argument("--ideal-method", default="PF_Durham_IdealMatching",
                        help="method dir for the ideal-matching (+ Jet definition) curve")
    parser.add_argument("--reco-hist", default="h_mH_reco",
                        help="jet-based mH histogram name (use h_mH_reco_fixed for a re-run "
                             "with the corrected H-jet mapping)")
    parser.add_argument("--rebin", type=int, default=1,
                        help="merge N adjacent 0.25 GeV bins for drawing (stats stay full-resolution)")
    parser.add_argument("--zoom-range", nargs=2, type=float, default=[90.0, 150.0],
                        metavar=("MIN", "MAX"))
    parser.add_argument("--wide-range", nargs=2, type=float, default=[60.0, 200.0],
                        metavar=("MIN", "MAX"))
    parser.add_argument("--peak-range", nargs=2, type=float, default=[105.0, 145.0],
                        metavar=("MIN", "MAX"),
                        help="tight window (mH +/- 20 GeV) used for "
                             "mH_jet_multiplicity_zoom.pdf (default: 105 145)")
    parser.add_argument("--processes", nargs="*", default=None,
                        help="subset of processes (default: all in process_config)")
    parser.add_argument("--compare-fixed", action="store_true",
                        help="additionally overlay h_mH_reco vs h_mH_reco_fixed per process, "
                             "sizing up the parton->reco-jet mapping bug on identical events "
                             "(needs a histmaker run that produced h_mH_reco_fixed)")
    args = parser.parse_args()

    input_dir = args.inputDir or default_input_dir()
    args.inputDir = input_dir
    output_dir = args.outputDir or os.path.join(input_dir, "plots", "mh_decomposition")
    os.makedirs(output_dir, exist_ok=True)
    method_dirs = {"reco": args.reco_method, "ideal": args.ideal_method}
    processes = args.processes or list(PROCESS_TO_ROW_COL)

    print(f"input:  {input_dir}")
    print(f"methods: reco={args.reco_method}  ideal={args.ideal_method}  hist={args.reco_hist}")
    print(f"output: {output_dir}")

    data, warnings = collect(input_dir, method_dirs, args.reco_hist, processes, args.rebin)
    for w in warnings:
        print(f"  !! {w}")
    if not data:
        raise SystemExit("no histograms could be read - check --inputDir / --reco-method")

    figures = []
    # The PFlow-only label schemes go in their own subfolder so they don't get
    # mixed up with the gen-momentum ladders.
    pflow_dir = os.path.join(output_dir, "pflow")
    for variant_list, ladders, dest in ((VARIANTS, LADDERS, output_dir),
                                        (VARIANTS_SIMPLE, SIMPLE_LADDERS, output_dir),
                                        (VARIANTS_GENPHYS, GENPHYS_LADDERS, output_dir),
                                        (VARIANTS_PFLOW, PFLOW_LADDERS, pflow_dir)):
        os.makedirs(dest, exist_ok=True)
        for stem, keys in ladders.items():
            n = len(keys)
            for suffix, x_range, log_y in (
                ("", args.zoom_range, False),
                ("_wide_log", args.wide_range, True),
            ):
                for t_suffix, transposed in (("", False), ("_by_jets", True)):
                    name = f"{stem}{t_suffix}{suffix}.pdf"
                    fig = fig_grid(data, keys, method_dirs, args.reco_hist,
                                   tuple(x_range), log_y, all_variants=variant_list,
                                   transposed=transposed)
                    fig.savefig(os.path.join(dest, name))
                    plt.close(fig)
                    scale = "log y" if log_y else "linear"
                    orient = ("2/4/6 jets across the row, flavour down the rows"
                              if transposed else
                              "flavour across the row, 2/4/6 jets down the rows")
                    rel = os.path.relpath(os.path.join(dest, name), output_dir)
                    figures.append((rel, f"{n} curves, {x_range[0]:g}-{x_range[1]:g} GeV, "
                                         f"{scale}; {orient}"))
                    print(f"  wrote {os.path.join(dest, name)}")

    # Poster version: genphys with the gluon row dropped, no guide arrows, tight.
    no_gluon_rows = [JETS_IN_COLUMNS[0], JETS_IN_COLUMNS[2]]
    no_gluon_processes = [pr for row in no_gluon_rows for pr in row]
    gluon_row = [JETS_IN_COLUMNS[1]]
    for poster_name, kwargs, orient in (
        ("mH_decomposition_genphys_by_jets_NoGluons.pdf",
         dict(transposed=True, jets_rows=no_gluon_rows),
         "2/4/6 jets across the row, flavour down the rows"),
        # rows <-> columns: 2/4/6 jets down the rows, flavour across
        ("mH_decomposition_genphys_by_jets_NoGluons_T.pdf",
         dict(transposed=False, only_processes=no_gluon_processes),
         "2/4/6 jets down the rows, flavour across the row"),
        # the gluon row on its own, the counterpart of the NoGluons figure
        ("mH_decomposition_genphys_by_jets_GluonsOnly.pdf",
         dict(transposed=True, jets_rows=gluon_row),
         "gluon final states only, 2/4/6 jets across the row"),
    ):
        fig = fig_grid(data, GENPHYS_LADDERS["mH_decomposition_genphys"], method_dirs,
                       args.reco_hist, tuple(args.zoom_range), False,
                       all_variants=VARIANTS_GENPHYS, arrows=False,
                       jet_badge=True, **kwargs)
        fig.savefig(os.path.join(output_dir, poster_name))
        plt.close(fig)
        figures.append((poster_name,
                        "poster version: genphys, light-flavour and b-rich only "
                        f"(no gluon row), {orient}, no guide arrows"))
        print(f"  wrote {os.path.join(output_dir, poster_name)}")

    zoom_view = (tuple(args.zoom_range), False)
    wide_view = (tuple(args.wide_range), True)
    peak_view = (tuple(args.peak_range), False)

    # Poster/talk figures: one process per page, the full ladder, big fonts.
    pp_path, pp_pages = write_per_process_pdf(
        data, output_dir, LADDERS["mH_decomposition_all"], [peak_view, wide_view])
    figures.append((os.path.basename(pp_path),
                    f"{pp_pages} pages: one poster-sized canvas per process (ordered "
                    f"2-jet -> 4-jet -> 6-jet) with the full ladder, peak-zoom "
                    f"({args.peak_range[0]:g}-{args.peak_range[1]:g} GeV) and log-y wide page each"))
    print(f"  wrote {pp_path} ({pp_pages} pages)")

    for name, views, desc in (
        ("mH_jet_multiplicity.pdf", [zoom_view, wide_view],
         "linear zoom and log-y wide page for each"),
        ("mH_jet_multiplicity_zoom.pdf", [(tuple(args.peak_range), False)],
         f"linear, tight peak window {args.peak_range[0]:g}-{args.peak_range[1]:g} GeV"),
    ):
        mult_path, mult_pages = write_multiplicity_pdf(data, output_dir, views, name=name)
        figures.append((name,
                        f"{mult_pages} pages: one canvas per mH definition with the 2/4/6-jet "
                        f"processes overlaid, light-flavour and b-jet-rich sets; {desc}"))
        print(f"  wrote {mult_path} ({mult_pages} pages)")

    rows = table_rows(data)

    # Common-selection version: only produced if the histmaker run that made
    # these histograms had both jet definitions (h_mH_common_* present).
    for cvars, cladders in ((PFLOW_COMMON_VARIANTS, PFLOW_COMMON_LADDERS),
                            (GENPHYS_COMMON_VARIANTS, GENPHYS_COMMON_LADDERS)):
      common, common_warn = collect(input_dir, method_dirs, args.reco_hist, processes,
                                  args.rebin, variants=cvars)
      if common:
        keys = [v["key"] for v in cvars]
        for stem, kk in cladders.items():
            for suffix, x_range, log_y in (("", args.zoom_range, False),
                                           ("_wide_log", args.wide_range, True)):
                for t_suffix, transposed in (("", False), ("_by_jets", True)):
                    name = f"{stem}{t_suffix}{suffix}.pdf"
                    fig = fig_grid(common, kk, method_dirs, args.reco_hist,
                                   tuple(x_range), log_y,
                                   all_variants=cvars,
                                   transposed=transposed)
                    os.makedirs(pflow_dir, exist_ok=True)
                    fig.savefig(os.path.join(pflow_dir, name))
                    plt.close(fig)
                    figures.append((os.path.join("pflow", name),
                                    "PFlow definitions on a common event "
                                    "sample (Higgs jets found by both jet "
                                    "definitions)"))
                    print(f"  wrote {os.path.join(pflow_dir, name)}")
        rows += table_rows(common, variants=cvars)
      else:
        print("  (no h_mH_common_* histograms in this tree - skipping the "
              "common-selection figures)")

    for path in write_tables(rows, output_dir):
        print(f"  wrote {path}")

    if args.compare_fixed:
        figures += fix_comparison(input_dir, method_dirs, args, processes, output_dir)

    write_readme(output_dir, args, method_dirs, figures, data)
    print(f"  wrote {os.path.join(output_dir, 'README.md')}")


if __name__ == "__main__":
    main()
