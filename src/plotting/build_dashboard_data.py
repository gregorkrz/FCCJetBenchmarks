"""Stage 3 of the resolution pipeline (no ROOT, no fitting — pure packaging).

Walks a top-level histogram directory containing one subfolder per jet
clustering method (e.g. `PF_Durham`, `CaloJets_Durham`, `PF_AntiKtR08`, ...),
loads each method's `plots_resolution/resolution_dashboard_data.pkl`
(produced by `resolution_plots.py`) and `plots_mass/mass_dashboard_data.pkl`
(produced by `mass_plots.py`), and consolidates everything for
`make_interactive_dashboard.py` into:

- `dashboard_data.json`: the light, downsampled data that gets inlined into
  dashboard.html itself. Every full-resolution histogram
  (`_extract_full_histograms`) is replaced by a `"full_url"` pointer.
- `full_hist/...json`: one small JSON file per full-resolution histogram
  (one per energy/angle/eta/costheta bin, one per mass definition), so
  dashboard.html can fetch, on demand, only the exact histogram the user
  clicks on instead of the entire full-resolution dataset.

Usage:
    python src/plotting/build_dashboard_data.py --inputDir $PATH_TO_HISTOGRAMS
"""
import argparse
import json
import os
import pickle
import re
import shutil

import numpy as np

from src.process_config import (
    HUMAN_READABLE_PROCESS_NAMES,
    LINE_STYLES,
    NUMBER_OF_JETS,
    PROCESS_COLORS,
    PROCESS_TO_ROW_COL,
)
from src.plotting.resolution_methods import (
    RESOLUTION_MODELS,
    fit_resolution_model,
    fit_resolution_staged_const,
)

# Reproduce the energy-JER bounds used in resolution_plots.py so the weighted
# refit here stays consistent with the unweighted fit produced upstream.
ENERGY_JER_BOUNDS = ([0.0, 0.005, 0.0], [np.inf, 0.04, np.inf])

DASHBOARD_PKL_RELPATH = os.path.join("plots_resolution", "resolution_dashboard_data.pkl")
MASS_DASHBOARD_PKL_RELPATH = os.path.join("plots_mass", "mass_dashboard_data.pkl")
STATS_JSON_RELPATH = "basic_stats_summary.json"


def method_label(method_name):
    """Turn a method directory name into a human-readable label.

    Mirrors the naming conventions used in joint_plots.py's method_dict,
    generalized to handle arbitrary anti-kt radii.
    """
    if method_name == "PF_Durham":
        return "PF Durham"
    if method_name == "PF_Durham_IdealMatching":
        return "PF Durham (Ideal Matching)"
    if method_name == "CaloJets_Durham":
        return "Calo Durham"

    m = re.match(r"^PF_E_recovery_AntiKtR(\d+)$", method_name)
    if m:
        return f"PF AntiKt R={m.group(1)} (E recovery)"

    m = re.match(r"^PF_AntiKtR(\d+)$", method_name)
    if m:
        return f"PF AntiKt R={m.group(1)}"

    return method_name


def _sanitize(obj):
    """Recursively convert numpy types/arrays to plain JSON-serializable Python types."""
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _sanitize(obj.tolist())
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    return obj


def _sanitize_path_part(part):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(part))


def _extract_full_histograms(obj, path_parts, full_dir, url_prefix):
    """Pull the full-resolution "edges_full"/"y_full" pair out of every
    histogram-bin-like dict in obj, writing each one to its own small JSON
    file under full_dir (named after the path to reach it, e.g.
    `PF_Durham/p8_ee_ZH_qqqq_ecm240/energy/_all/bins/3.json`) and replacing it
    in the returned light copy with a `"full_url"` pointer to that file - so
    dashboard.html can lazily fetch exactly the one histogram it needs to
    render, instead of the entire full-resolution dataset.
    """
    if isinstance(obj, dict):
        if "edges_full" in obj and "y_full" in obj:
            light = {k: v for k, v in obj.items() if k not in ("edges_full", "y_full")}
            rel_path = "/".join(_sanitize_path_part(p) for p in path_parts) + ".json"
            full_path = os.path.join(full_dir, rel_path)
            os.makedirs(os.path.dirname(full_path), exist_ok=True)
            with open(full_path, "w") as fh:
                json.dump({"edges": obj["edges_full"], "y": obj["y_full"]}, fh)
            light["full_url"] = url_prefix + rel_path
            return light
        return {
            k: _extract_full_histograms(v, path_parts + [k], full_dir, url_prefix)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [
            _extract_full_histograms(v, path_parts + [i], full_dir, url_prefix)
            for i, v in enumerate(obj)
        ]
    return obj


def _point_errors(mid_points, values, bins):
    """Analytic per-point 1-sigma error on sigma/E for each energy bin.

    Option A: the statistical error on a width estimator from N events scales
    as delta(sigma)/sigma ~ 1/sqrt(2N) (the Gaussian-sigma result; the
    narrowest-68% interval shares the same 1/sqrt(N) scaling). The MPV
    (denominator of sigma/E) has a much smaller relative error, so to leading
    order delta(sigma/E) ~ (sigma/E)/sqrt(2N). N per point is the n_jets stored
    alongside each bin record. Returns None if the bins don't line up with the
    fitted points (older data) so we fall back to the unweighted fit.
    """
    if not bins or len(bins) != len(values):
        return None
    errors = []
    for r, b in zip(values, bins):
        n = b.get("n_jets", 0) if isinstance(b, dict) else 0
        errors.append(abs(r) / np.sqrt(2.0 * n) if n and n > 0 else np.nan)
    return errors


def _chi2_ndf(model, mid_points, values, errors, popt):
    """Reduced chi2 of a weighted fit over the points that carried a valid error."""
    func = RESOLUTION_MODELS[model]["func"]
    x = np.asarray(mid_points, dtype=float)
    y = np.asarray(values, dtype=float)
    e = np.asarray(errors, dtype=float)
    good = np.isfinite(e) & (e > 0)
    ndf = int(good.sum()) - len(popt)
    if ndf <= 0:
        return None
    resid = (y[good] - func(x[good], *popt)) / e[good]
    return float(np.sum(resid ** 2) / ndf)


# Energy-JER fit variants to add alongside the unweighted three_param/two_param
# produced upstream by resolution_plots.py. Keyed by the fit name shown in the
# dashboard -> (base model in RESOLUTION_MODELS, extra fit_resolution_model
# kwargs). This is the "top 10" set selected by the fit_trials/ experiments,
# ranked by leave-one-out CV MAE (see fit_trials/README.md), each offered in
# unweighted and stat-error-weighted (Option A) form so they sit side by side
# with the incumbents in the fit-model dropdown.
#
# Ordering roughly follows CV rank. `logE_4param` (weighted) was the overall
# CV winner (~3x better than the bounded production fit, and unlike the 5-param
# forms its CV error ~= its in-sample error, i.e. it generalises rather than
# overfits). `logE_and_invE2` has the best in-sample MAE but overfits under CV -
# included (both forms) so that trade-off is visible on the site.
WEIGHTED_FIT_VARIANTS = {
    # rank 1 (CV): extra d*log(E) drift term - the standout generaliser.
    "logE_4param_weighted": ("logE_4param", {}),
    "logE_4param": ("logE_4param", {"errors": None}),
    # extra steep low-E terms (great in-sample, solid CV).
    "invE15_noise": ("invE15_noise", {"errors": None}),
    "invE2_noise": ("invE2_noise", {"errors": None}),
    "invE2_noise_weighted": ("invE2_noise", {}),
    # no-stochastic-term 4-param (parsimony variant).
    "logE_invE2_noStoch": ("logE_invE2_noStoch", {"errors": None}),
    # 5-param kitchen sink (best in-sample; overfits under CV - shown for compare).
    "logE_and_invE2": ("logE_and_invE2", {"errors": None}),
    "logE_and_invE2_weighted": ("logE_and_invE2", {}),
    # the incumbent three_param, now also offered weighted (Option A).
    "three_param_weighted": ("three_param", {"bounds_override": ENERGY_JER_BOUNDS}),
    "two_param_weighted": ("two_param", {}),
    # the staged "pin B from the high-E tail, then fit A and C" strategy
    # (fit_trials/ trial 03). Same three_param algebra, different fitting
    # procedure, hence the custom `fitter`; it ranks last on CV but is offered so
    # the strategy can be inspected per method/process rather than only in the
    # aggregate leaderboard.
    "staged_const_first": (
        "three_param",
        {"errors": None, "fitter": fit_resolution_staged_const},
    ),
}


def _add_weighted_energy_fits(sanitized_method_data):
    """Add Option A weighted-fit variants to every process's energy entry.

    Mutates the per-process dict in place. For each energy resolution entry it
    computes delta(sigma/E) = (sigma/E)/sqrt(2N) from the stored per-bin n_jets,
    attaches it as sigma_over_E_err (for the dashboard error bars), and adds new
    entries to `fits` - "three_param_weighted"/"two_param_weighted" - fitted with
    those errors (sigma=, absolute_sigma=True), each carrying popt/popt_err/
    chi2_ndf. The existing unweighted fits and the entry's default popt/fit_x/
    fit_y are intentionally left as-is. Only the energy quantity is touched;
    angular fits are unchanged.
    """
    for process, proc_data in sanitized_method_data.items():
        energy = proc_data.get("energy") if isinstance(proc_data, dict) else None
        if not isinstance(energy, dict):
            continue
        for jet_part, entry in energy.items():
            if not isinstance(entry, dict):
                continue
            mids = entry.get("mid_points")
            vals = entry.get("sigma_over_E")
            if not mids or not vals or len(mids) < 2:
                continue
            errors = _point_errors(mids, vals, entry.get("bins"))
            if errors is not None:
                entry["sigma_over_E_err"] = errors
            fits = entry.setdefault("fits", {})
            for fit_name, (base_model, kw) in WEIGHTED_FIT_VARIANTS.items():
                kw = dict(kw)
                # A variant may set errors=None to force an unweighted fit;
                # otherwise it inherits the analytic per-point errors.
                fit_errors = kw.pop("errors", errors)
                # ...and may supply its own fitter (a staged/multi-step
                # procedure) as long as it returns fit_resolution_model's
                # (xs, ys, popt, pcov) with popt in `base_model` order.
                fitter = kw.pop("fitter", fit_resolution_model)
                try:
                    xs, ys, popt, pcov = fitter(
                        mids, vals, model=base_model, errors=fit_errors, **kw
                    )
                except Exception as e:
                    print(f"    {fit_name} fit failed for {process}/{jet_part}: {e}")
                    continue
                fit_out = {
                    "popt": popt.tolist(),
                    "fit_x": xs.tolist(),
                    "fit_y": ys.tolist(),
                    "popt_err": np.sqrt(np.abs(np.diag(pcov))).tolist(),
                    "base_model": base_model,
                }
                # chi2/ndf always uses the statistical errors (comparable across
                # weighted and unweighted variants), when available.
                if errors is not None:
                    fit_out["chi2_ndf"] = _chi2_ndf(base_model, mids, vals, errors, popt)
                fits[fit_name] = fit_out


def discover_methods(input_dir):
    methods = []
    for name in sorted(os.listdir(input_dir)):
        if name.startswith("plots"):
            continue
        method_dir = os.path.join(input_dir, name)
        pkl_path = os.path.join(method_dir, DASHBOARD_PKL_RELPATH)
        if os.path.isdir(method_dir) and os.path.isfile(pkl_path):
            methods.append((name, pkl_path))
    return methods


def load_mass_data(input_dir, method_name):
    """Load a method's Higgs-mass dashboard pickle, if it exists.

    Independent of discover_methods() so a method missing mass_plots.py
    output still shows up in the dashboard, just without a mass quantity.
    """
    pkl_path = os.path.join(input_dir, method_name, MASS_DASHBOARD_PKL_RELPATH)
    if not os.path.isfile(pkl_path):
        return None
    with open(pkl_path, "rb") as fh:
        return pickle.load(fh)


def load_stats(input_dir):
    """Load the basic_stats_summary.json produced by print_basic_stats.py, if present."""
    stats_path = os.path.join(input_dir, STATS_JSON_RELPATH)
    if not os.path.isfile(stats_path):
        return None
    with open(stats_path) as fh:
        return json.load(fh)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputDir", type=str, required=True)
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON path (default: $inputDir/plots/dashboard_data.json)",
    )
    parser.add_argument(
        "--full-hist-dir",
        type=str,
        default=None,
        help="Output directory for the per-histogram full-resolution JSON files, lazily fetched "
        "by dashboard.html (default: alongside --output, named full_hist/)",
    )
    args = parser.parse_args()

    methods = discover_methods(args.inputDir)
    if not methods:
        raise FileNotFoundError(
            f"No method subfolders with {DASHBOARD_PKL_RELPATH} found under {args.inputDir}. "
            "Run extract_resolution_data.py and resolution_plots.py for each method first."
        )

    output_path = args.output or os.path.join(args.inputDir, "plots", "dashboard_data.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    full_hist_dir = args.full_hist_dir or os.path.join(os.path.dirname(output_path), "full_hist")
    if os.path.isdir(full_hist_dir):
        shutil.rmtree(full_hist_dir)
    url_prefix = os.path.relpath(full_hist_dir, os.path.dirname(output_path)) + "/"

    all_processes = set()
    methods_out = {}
    for method_name, pkl_path in methods:
        print("Loading dashboard data for method:", method_name)
        with open(pkl_path, "rb") as fh:
            data = pickle.load(fh)
        all_processes.update(data.keys())

        mass_data = load_mass_data(args.inputDir, method_name)
        if mass_data is not None:
            print("  + Higgs mass dashboard data for method:", method_name)
            for process, mass_entry in mass_data.items():
                data.setdefault(process, {})["mass"] = mass_entry

        sanitized = _sanitize(data)
        # Option A: add weighted energy-JER fit variants (plus per-point errors /
        # param uncertainties / chi2) alongside the existing unweighted fits,
        # before the full histograms are split out - the bin records still carry
        # n_jets at this point.
        _add_weighted_energy_fits(sanitized)
        light_data = _extract_full_histograms(
            sanitized, [method_name], full_hist_dir, url_prefix
        )
        methods_out[method_name] = {
            "label": method_label(method_name),
            "processes": light_data,
        }

    process_meta = {}
    for process in sorted(all_processes):
        # PROCESS_TO_ROW_COL stores (col, row): col = jet-multiplicity band
        # (2/4/6 jets, "higher number of final-state jets" axis), row = B-hadron
        # content band. Used by the dashboard's Grid tab to lay processes out in
        # the same 5x3 matrix as joint_plots.py.
        grid_col, grid_row = PROCESS_TO_ROW_COL.get(process, (None, None))
        process_meta[process] = {
            "label": HUMAN_READABLE_PROCESS_NAMES.get(process, process),
            "color": PROCESS_COLORS.get(process),
            "line_style": LINE_STYLES.get(process, "-"),
            "n_jets": NUMBER_OF_JETS.get(process),
            "grid_row": grid_row,
            "grid_col": grid_col,
        }

    stats = load_stats(args.inputDir)
    if stats is not None:
        print("Loaded basic_stats_summary.json")
    else:
        print(f"No {STATS_JSON_RELPATH} found under {args.inputDir}; Statistics tab will be empty.")

    output = {
        "methods": methods_out,
        "process_meta": process_meta,
        "stats": stats or {},
    }

    with open(output_path, "w") as fh:
        json.dump(output, fh)

    print("Saved consolidated dashboard data to:", output_path)
    print(f"Methods: {len(methods_out)}, Processes: {len(process_meta)}")

    if os.path.isdir(full_hist_dir):
        n_full_histograms = sum(len(files) for _, _, files in os.walk(full_hist_dir))
        print(
            f"Saved {n_full_histograms} full-resolution histograms to:", full_hist_dir,
            "(one JSON file per histogram, lazily fetched by dashboard.html only for the exact "
            "histogram being displayed when the 'show full histogram' checkbox is toggled - keep "
            "it next to dashboard.html and serve both over http(s), not file://)",
        )


if __name__ == "__main__":
    main()
