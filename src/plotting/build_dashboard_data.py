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
)

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

        light_data = _extract_full_histograms(
            _sanitize(data), [method_name], full_hist_dir, url_prefix
        )
        methods_out[method_name] = {
            "label": method_label(method_name),
            "processes": light_data,
        }

    process_meta = {}
    for process in sorted(all_processes):
        process_meta[process] = {
            "label": HUMAN_READABLE_PROCESS_NAMES.get(process, process),
            "color": PROCESS_COLORS.get(process),
            "line_style": LINE_STYLES.get(process, "-"),
            "n_jets": NUMBER_OF_JETS.get(process),
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
