"""Higgs-mass grid figures for the three generalized e+e- exponents.

Two complementary figure sets over the same 5x3 process matrix that
joint_plots.py uses (`PROCESS_TO_ROW_COL`):

  per-exponent   one figure per exponent, overlaying that exponent's six radii
                 plus Durham. Algorithm held fixed, radius varied.
                 -> mH_grid_exponent_p{m1,0,p1}.pdf

  per-radius     one figure per radius, overlaying p = -1, 0, +1 plus Durham.
                 Radius held fixed, algorithm varied. This is the direct
                 algorithm comparison and is not expressible through
                 joint_plots.py's one-family-at-a-time `--family` mode.
                 -> mH_grid_radius_R{04..14}.pdf

Unlike joint_plots.py this reads *only* `plots_mass/Higgs_mass_histograms_data.pkl`,
so it needs `mass_plots.py` to have run but not the two much slower resolution
stages (`extract_resolution_data.py` + `resolution_plots.py`).

Usage:
    source env.sh
    python src/plotting/mh_grid_plots.py --inputDir $PATH_TO_HISTOGRAMS
"""
import argparse
import os
import pickle
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from src.process_config import (  # noqa: E402
    HUMAN_READABLE_PROCESS_NAMES,
    NUMBER_OF_JETS,
    PROCESS_TO_ROW_COL,
    RADIUS_SCAN,
    family_prefixes,
    radius_to_str,
)

MASS_PKL_RELPATH = os.path.join("plots_mass", "Higgs_mass_histograms_data.pkl")
REFERENCE_METHOD = "PF_Durham"
REFERENCE_LABEL = "Durham"
REFERENCE_COLOR = "#1f77b4"

# The three exponents of JetClustering::clustering_ee_genkt, in the order they
# should appear in a legend. `family` keys match src/process_config.JET_FAMILIES;
# the colours match the convention used across the repo's figures and slides.
EXPONENTS = [
    dict(p=-1, family="ee-akt", tag="m1", label=r"ee-anti-$k_T$",
         color="#2A9D8F", cmap=plt.cm.GnBu, cmap_range=(0.40, 0.92)),
    dict(p=0, family="ee-ca", tag="0", label="ee-C/A",
         color="#7A3CBF", cmap=plt.cm.Purples, cmap_range=(0.40, 0.90)),
    dict(p=1, family="ee-kt", tag="p1", label=r"ee-$k_T$",
         color="#D18B00", cmap=plt.cm.YlOrBr, cmap_range=(0.35, 0.90)),
]

# With inclusive clustering a 6-jet event merges into fewer jets as R grows, so
# the fully-matched filter leaves almost nothing at large R (measured for C/A:
# ~1.4% at R=1.2, ~0.001% at R=1.4) and the surviving mH curve is pure noise.
# Curves above this radius are dropped from the 6-jet panels; every drop is
# reported on stdout rather than being silent.
MAX_RADIUS_FOR_6JETS = 1.0
MH_XLIM = (90, 150)


def fmt_p(p):
    """LaTeX-ready exponent label: -1, 0, +1 (not "+0")."""
    return f"{p:+d}" if p != 0 else "0"


def rgba_to_hex(rgba):
    r, g, b = rgba[:3]
    return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))


def resolve_method(input_dir, family, radius):
    """Directory for one (family, radius), trying canonical then legacy names."""
    for prefix in family_prefixes(family):
        candidate = f"{prefix}{radius_to_str(radius)}"
        if os.path.isfile(os.path.join(input_dir, candidate, MASS_PKL_RELPATH)):
            return candidate
    return None


def load_mass(input_dir, method):
    with open(os.path.join(input_dir, method, MASS_PKL_RELPATH), "rb") as fh:
        return pickle.load(fh)


def pass_rate(input_dir, method, process):
    """Fraction of events surviving the fully-matched-jets filter, or None.

    Written per (method, process) by histmaker.py as basic_stats_<process>.pkl.
    Shown in the legend because it varies enormously across the radius scan -
    at large R an inclusive 6-jet event merges into fewer jets and the filter
    keeps ~1% of events, so a curve can look reasonable while resting on almost
    no statistics.
    """
    path = os.path.join(input_dir, method, f"basic_stats_{process}.pkl")
    try:
        with open(path, "rb") as fh:
            stats = pickle.load(fh)
        before = stats["before_filtering"]
        return stats["after_filtering"] / before if before > 0 else None
    except Exception:
        return None


def step_hist(ax, entry, label, color, linestyle):
    """Draw one mH histogram as a step outline. Returns the integral."""
    x = np.asarray(entry["x_vals_reco"], dtype=float)
    y = np.asarray(entry["y_vals_reco"], dtype=float)
    if x.size < 2:
        return 0.0
    width = x[1] - x[0]
    edges = np.concatenate((x - width / 2, [x[-1] + width / 2]))
    ax.hist(x, bins=edges, weights=y, histtype="step",
            label=label, color=color, linestyle=linestyle)
    return float(y.sum())


def new_grid():
    fig, ax = plt.subplots(5, 3, figsize=(10, 11))
    return fig, ax


def finish_grid(fig, ax, title):
    for row in range(ax.shape[0]):
        for col in range(ax.shape[1]):
            a = ax[row, col]
            if not a.has_data():
                a.set_axis_off()
                continue
            a.set_xlabel("$m_H$ [GeV]")
            a.set_xlim(*MH_XLIM)
            a.grid()
            a.legend(fontsize=6.0, loc="upper left", framealpha=0.85,
                     borderpad=0.3, labelspacing=0.25, handlelength=1.6)
    fig.suptitle(title + "\n" + r"legend: (n) = fraction of events passing the "
                 r"fully-matched-jets filter for that process",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.975))


def draw_panels(ax, curves, input_dir):
    """curves: list of (method, label, color, linestyle). Returns dropped list."""
    dropped = []
    caches = {}
    for method, label, color, linestyle in curves:
        try:
            caches[method] = load_mass(input_dir, method)
        except Exception as exc:
            dropped.append(f"{method} (unreadable: {type(exc).__name__})")
    for process, (col, row) in PROCESS_TO_ROW_COL.items():
        a = ax[row, col]
        drawn = False
        for method, label, color, linestyle in curves:
            data = caches.get(method)
            if not data or process not in data:
                continue
            radius = None
            for e in EXPONENTS:
                r = _radius_of(method, e["family"])
                if r is not None:
                    radius = r
                    break
            if (radius is not None and NUMBER_OF_JETS.get(process) == 6
                    and radius > MAX_RADIUS_FOR_6JETS):
                dropped.append(f"{method} from {process} (R={radius} > "
                               f"{MAX_RADIUS_FOR_6JETS} in a 6-jet panel)")
                continue
            # Legend carries this method's filter pass rate *for this process*,
            # so a curve resting on 1% of the events is visibly flagged.
            rate = pass_rate(input_dir, method, process)
            panel_label = label if rate is None else f"{label} ({rate:.2f})"
            step_hist(a, data[process], panel_label, color, linestyle)
            drawn = True
        if not drawn:
            continue
        # Jet count in the title rather than as a floating badge: with five
        # legend entries per panel an in-axes badge collides with the legend.
        a.set_title(
            "{}   [{} jets]".format(
                HUMAN_READABLE_PROCESS_NAMES.get(process, process),
                NUMBER_OF_JETS[process]),
            fontsize=9)
    return dropped


def _radius_of(method, family):
    for prefix in family_prefixes(family):
        if method.startswith(prefix):
            digits = method[len(prefix):]
            if digits.isdigit():
                return int(digits) / 10
    return None


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputDir", default=os.environ.get("PATH_TO_HISTOGRAMS"),
                   help="Histogram tree root (default $PATH_TO_HISTOGRAMS)")
    p.add_argument("--outputDir", default=None,
                   help="Default: <inputDir>/plots/mh_grids")
    p.add_argument("--sets", default="exponent,radius",
                   help="Which figure sets to draw: exponent, radius, or both")
    args = p.parse_args(argv)
    if not args.inputDir:
        raise SystemExit("set PATH_TO_HISTOGRAMS or pass --inputDir")
    input_dir = args.inputDir.rstrip("/")
    output_dir = args.outputDir or os.path.join(input_dir, "plots", "mh_grids")
    wanted = {s.strip() for s in args.sets.split(",") if s.strip()}

    have_reference = os.path.isfile(
        os.path.join(input_dir, REFERENCE_METHOD, MASS_PKL_RELPATH))
    if not have_reference:
        print(f"note: {REFERENCE_METHOD} has no {MASS_PKL_RELPATH}; "
              "the Durham reference curve will be omitted")

    # Which (exponent, radius) combinations actually have a mass pickle.
    available = {}
    for e in EXPONENTS:
        for radius in RADIUS_SCAN:
            method = resolve_method(input_dir, e["family"], radius)
            if method:
                available[(e["p"], radius)] = method
    if not available:
        raise SystemExit(
            f"no radius-scan method directory under {input_dir} has "
            f"{MASS_PKL_RELPATH} yet - run mass_plots.py on them first")

    print(f"found {len(available)} of {len(EXPONENTS) * len(RADIUS_SCAN)} "
          "(exponent, radius) combinations with mass pickles:")
    for e in EXPONENTS:
        radii = [r for (pp, r) in sorted(available) if pp == e["p"]]
        print(f"  p={fmt_p(e['p']):>2s} {e['label']:16s} R = "
              + (", ".join(f"{r:.1f}" for r in radii) if radii else "(none)"))

    os.makedirs(output_dir, exist_ok=True)
    written, all_dropped = [], []

    # --- one figure per exponent: that exponent's radii, light -> dark ---------
    if "exponent" in wanted:
        for e in EXPONENTS:
            radii = [r for (pp, r) in sorted(available) if pp == e["p"]]
            if not radii:
                print(f"skipping per-exponent grid for p={fmt_p(e['p'])}: no radii present")
                continue
            shades = np.linspace(*e["cmap_range"], max(1, len(radii)))
            curves = []
            if have_reference:
                curves.append((REFERENCE_METHOD, REFERENCE_LABEL, REFERENCE_COLOR, "-"))
            for i, radius in enumerate(radii):
                curves.append((available[(e["p"], radius)],
                               f"{e['label']} R={radius:.1f}",
                               rgba_to_hex(e["cmap"](shades[i])), "-"))
            fig, ax = new_grid()
            all_dropped += draw_panels(ax, curves, input_dir)
            finish_grid(fig, ax,
                        f"$m_H$ per process - {e['label']} radius scan "
                        f"($p = {fmt_p(e['p'])}$), vs Durham")
            path = os.path.join(output_dir, f"mH_grid_exponent_p{e['tag']}.pdf")
            fig.savefig(path)
            plt.close(fig)
            written.append(path)

    # --- one figure per radius: the three exponents at that radius ------------
    if "radius" in wanted:
        for radius in RADIUS_SCAN:
            present = [e for e in EXPONENTS if (e["p"], radius) in available]
            if not present:
                print(f"skipping per-radius grid for R={radius:.1f}: no exponents present")
                continue
            curves = []
            if have_reference:
                curves.append((REFERENCE_METHOD, REFERENCE_LABEL, REFERENCE_COLOR, "-"))
            for e in present:
                curves.append((available[(e["p"], radius)],
                               f"{e['label']} ($p={fmt_p(e['p'])}$)", e["color"], "-"))
            fig, ax = new_grid()
            all_dropped += draw_panels(ax, curves, input_dir)
            finish_grid(fig, ax,
                        f"$m_H$ per process - the three exponents at "
                        f"$R = {radius:.1f}$, vs Durham")
            path = os.path.join(output_dir, f"mH_grid_radius_R{radius_to_str(radius)}.pdf")
            fig.savefig(path)
            plt.close(fig)
            written.append(path)

    if all_dropped:
        print(f"\ndropped {len(all_dropped)} curve(s):")
        for d in sorted(set(all_dropped)):
            print(f"  {d}")
    print(f"\nwrote {len(written)} figure(s) to {output_dir}:")
    for path in written:
        print("  " + os.path.basename(path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
