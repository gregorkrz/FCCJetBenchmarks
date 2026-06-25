"""Stage 3 of the resolution pipeline (no ROOT, no fitting — pure packaging).

Walks a top-level histogram directory containing one subfolder per jet
clustering method (e.g. `PF_Durham`, `CaloJets_Durham`, `PF_AntiKtR08`, ...),
loads each method's `plots_resolution/resolution_dashboard_data.pkl`
(produced by `resolution_plots.py`), and consolidates everything into a
single JSON file consumed by `make_interactive_dashboard.py`.

Usage:
    python src/plotting/build_dashboard_data.py --inputDir $PATH_TO_HISTOGRAMS
"""
import argparse
import json
import os
import pickle
import re

import numpy as np

from src.process_config import (
    HUMAN_READABLE_PROCESS_NAMES,
    LINE_STYLES,
    PROCESS_COLORS,
)

DASHBOARD_PKL_RELPATH = os.path.join("plots_resolution", "resolution_dashboard_data.pkl")


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputDir", type=str, required=True)
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON path (default: $inputDir/plots/dashboard_data.json)",
    )
    args = parser.parse_args()

    methods = discover_methods(args.inputDir)
    if not methods:
        raise FileNotFoundError(
            f"No method subfolders with {DASHBOARD_PKL_RELPATH} found under {args.inputDir}. "
            "Run extract_resolution_data.py and resolution_plots.py for each method first."
        )

    all_processes = set()
    methods_out = {}
    for method_name, pkl_path in methods:
        print("Loading dashboard data for method:", method_name)
        with open(pkl_path, "rb") as fh:
            data = pickle.load(fh)
        all_processes.update(data.keys())
        methods_out[method_name] = {
            "label": method_label(method_name),
            "processes": _sanitize(data),
        }

    process_meta = {}
    for process in sorted(all_processes):
        process_meta[process] = {
            "label": HUMAN_READABLE_PROCESS_NAMES.get(process, process),
            "color": PROCESS_COLORS.get(process),
            "line_style": LINE_STYLES.get(process, "-"),
        }

    output = {
        "methods": methods_out,
        "process_meta": process_meta,
    }

    output_path = args.output or os.path.join(args.inputDir, "plots", "dashboard_data.json")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(output, fh)

    print("Saved consolidated dashboard data to:", output_path)
    print(f"Methods: {len(methods_out)}, Processes: {len(process_meta)}")


if __name__ == "__main__":
    main()
