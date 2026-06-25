"""Stage 1 of the resolution pipeline (ROOT-dependent).

Reads the per-bin histograms produced by histmaker.py (one ROOT file per
process in --inputDir) and dumps their raw content (bin edges + counts) to
a single pickle file, `plots_resolution/resolution_histograms.pkl`.

This is the only step in the resolution pipeline that needs ROOT / the
histogram ROOT files. Everything downstream (computing sigma/MPV per bin,
fitting the energy dependence, plotting, and building the interactive
dashboard) only needs that pickle, runs without ROOT, and is cheap enough
to run locally.

Usage:
    python src/plotting/extract_resolution_data.py --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
"""
import argparse
import os
import pickle

import numpy as np
import ROOT

from src.histmaker_tools.utils import neg_format

binsE = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
bins_eta = [-5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 5]
bins_costheta = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]

JET_PART_TO_HISTOGRAM_PREFIX = {
    "_charged": "binned_E_Charged_reco_over_true_FullGenJet_",
    "_neutral": "binned_E_Neutral_reco_over_true_FullGenJet_",
    "_all": "binned_E_reco_over_true_",
    "_photons": "binned_E_Photon_reco_over_true_FullGenJet_",
}

ANGLE_TO_HISTOGRAM_PREFIX = {
    "theta": "binned_deltaTheta_",
    "phi": "binned_deltaPhi_",
    "eta": "binned_deltaEta_",
}


def root_file_get_hist_and_edges(root_file, hist_name):
    h = root_file.Get(hist_name)
    if not h:
        print(f"Warning: histogram {hist_name} not found")
        return None, None
    nb = h.GetNbinsX()
    edges = np.array(
        [h.GetXaxis().GetBinLowEdge(1)]
        + [h.GetXaxis().GetBinUpEdge(i) for i in range(1, nb + 1)]
    )
    y = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    assert len(edges) == len(y) + 1
    return y, edges


def extract_bins(root_file, hist_prefix, bins, suffix=""):
    """Extract the raw (edges, y) histogram for each [bins[i], bins[i+1]) bin."""
    out = []
    for i in range(len(bins) - 1):
        hist_name = f"{hist_prefix}{suffix}{neg_format(bins[i])}_{neg_format(bins[i + 1])}"
        y, edges = root_file_get_hist_and_edges(root_file, hist_name)
        out.append(
            {
                "lo": bins[i],
                "hi": bins[i + 1],
                "hist_name": hist_name,
                "edges": edges,
                "y": y,
            }
        )
    return out


def extract_process(root_path, jet_parts, do_angles=True):
    f = ROOT.TFile.Open(root_path)
    data = {"energy": {}, "angles": {}}
    for jet_part in jet_parts:
        data["energy"][jet_part] = extract_bins(
            f, JET_PART_TO_HISTOGRAM_PREFIX[jet_part], binsE
        )
    if do_angles:
        for angle, prefix in ANGLE_TO_HISTOGRAM_PREFIX.items():
            data["angles"][angle] = extract_bins(f, prefix, binsE)
    data["eta_binned_E"] = extract_bins(
        f, "binned_E_reco_over_true_", bins_eta, suffix="eta_"
    )
    data["costheta_binned_E"] = extract_bins(
        f, "binned_E_reco_over_true_", bins_costheta, suffix="costheta_"
    )
    f.Close()
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputDir", type=str, required=True)
    parser.add_argument(
        "--angles-only",
        action="store_true",
        help="Only extract the histograms needed for angular resolutions (skip per-jet-part energy histograms)",
    )
    args = parser.parse_args()

    input_dir = args.inputDir
    output_dir = os.path.join(input_dir, "plots_resolution")
    os.makedirs(output_dir, exist_ok=True)

    processes = sorted(
        f.replace(".root", "") for f in os.listdir(input_dir) if f.endswith(".root")
    )
    jet_parts = ["_all"] if args.angles_only else ["_photons", "_neutral", "_charged", "_all"]

    data = {}
    for process in processes:
        print("Extracting histograms for process:", process)
        root_path = os.path.join(input_dir, f"{process}.root")
        data[process] = extract_process(root_path, jet_parts, do_angles=True)

    out_path = os.path.join(output_dir, "resolution_histograms.pkl")
    with open(out_path, "wb") as fh:
        pickle.dump(
            {
                "processes": processes,
                "binsE": binsE,
                "bins_eta": bins_eta,
                "bins_costheta": bins_costheta,
                "jet_parts": jet_parts,
                "data": data,
            },
            fh,
        )
    print("Saved raw resolution histograms to:", out_path)


if __name__ == "__main__":
    main()
