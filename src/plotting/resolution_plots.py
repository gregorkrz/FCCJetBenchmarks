"""Stage 2 of the resolution pipeline (no ROOT dependency).

Loads the raw per-bin histograms dumped by extract_resolution_data.py
(`plots_resolution/resolution_histograms.pkl`), computes the resolution
(sigma_E/E, angular sigma, response) per energy/angle bin, fits the
energy/angle dependence, and produces all the resolution plots.

This step is cheap (numpy/scipy/matplotlib only) and is meant to be run
locally, including while iterating on fitting algorithms - see
resolution_methods.py for the pluggable SIGMA_METHODS / RESOLUTION_MODELS
registries used here.

In addition to the PDF plots, this script writes:
  - energy_fit_params_per_process.pkl / angle_fit_params_per_process.pkl
    (consumed by joint_plots.py, format unchanged)
  - resolution_dashboard_data.pkl: a richer, self-contained dump (fit
    points + per-bin histograms + fit curves) used by
    build_dashboard_data.py to build the interactive HTML dashboard.

Usage:
    python src/plotting/resolution_plots.py --inputDir $PATH_TO_HISTOGRAMS/METHOD_NAME
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from src.process_config import PROCESS_COLORS, HUMAN_READABLE_PROCESS_NAMES, LINE_STYLES
import pickle
import matplotlib
from src.plotting.resolution_methods import SIGMA_METHODS, fit_resolution_model

matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 12,
    }
)

parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", type=str, required=True)
parser.add_argument("--angles-only", action="store_true", help="Only plot angular resolutions")

args = parser.parse_args()

inputDir = args.inputDir

outputDir = os.path.join(args.inputDir, "plots_resolution")
os.makedirs(outputDir, exist_ok=True)
print("Saving to directory:", outputDir)

raw_data_path = os.path.join(outputDir, "resolution_histograms.pkl")
if not os.path.exists(raw_data_path):
    raise FileNotFoundError(
        f"{raw_data_path} not found. Run extract_resolution_data.py first:\n"
        f"  python src/plotting/extract_resolution_data.py --inputDir {inputDir}"
    )
with open(raw_data_path, "rb") as fh:
    raw = pickle.load(fh)

binsE = raw["binsE"]
bins_eta = raw["bins_eta"]
bins_costheta = raw["bins_costheta"]
processList = {p: {"fraction": 1} for p in raw["processes"]}


def print_params(popt):
    if len(popt) == 2:
        return f"A={round(popt[0], 2)} C={round(popt[1], 2)}"
    return f"A={round(popt[0], 2)} B={round(popt[2], 2)} C={round(popt[1], 2)}"


def downsample_for_dashboard(y, edges, max_points=150):
    """Rebin (y, edges) for storage in the (small) dashboard data file.

    Purely for visual display in the dashboard - not used for any
    statistical computation, which always happens on the full-resolution
    histogram loaded from resolution_histograms.pkl.
    """
    n = len(y)
    if n <= max_points:
        return y, edges
    factor = int(np.ceil(n / max_points))
    n_new = n // factor
    y_ds = np.array([y[i * factor : (i + 1) * factor].sum() for i in range(n_new)])
    edges_ds = np.array([edges[i * factor] for i in range(n_new)] + [edges[n_new * factor]])
    return y_ds, edges_ds


def compute_resolution_for_process(
    bins_data,
    sigma_method="std68",
    wmin=0.7,
    wmax=1.2,
    divide_by_MPV=True,
    x_label=r"$\frac{E_{reco}}{E_{true}}$",
    n_jets_in_bin_min=70000,
):
    sigma_func = SIGMA_METHODS[sigma_method]
    fig_hist, ax_hist = plt.subplots(3, 1, figsize=(8, 8.5))

    bin_mid_points = []
    lo_hi_MPV = []
    sigmaEoverE = []
    responses = []
    bins_to_histograms = {}
    bin_records = []
    field_names = []
    total_statistics = 0

    for i, b in enumerate(bins_data):
        field_names.append(b["hist_name"])
        y, edges = b["y"], b["edges"]
        if y is None:
            print(f"Skipping bin [{b['lo']}, {b['hi']}] due to missing histogram")
            continue

        bin_widths = np.diff(edges)
        area = np.sum(y * bin_widths)
        n_jets_in_bin = int(np.sum(y))
        total_statistics += n_jets_in_bin
        y_normalized = y / area if area != 0 else y

        if n_jets_in_bin > 50000:
            for ax in ax_hist:
                ax.step(
                    edges[:-1],
                    y_normalized,
                    where="post",
                    label=f"[{b['lo']}, {b['hi']}] GeV (N={n_jets_in_bin})",
                )
        bins_to_histograms[i] = [y_normalized, edges]

        result = sigma_func(y, edges, wmin=wmin, wmax=wmax)
        if result is None:
            print(f"⚠️ {sigma_method} failed for bin [{b['lo']}, {b['hi']}]; aborting process.")
            return [], [], fig_hist, [], {}, [], [], 0, []
        std68, low, high, MPV = result

        lo_hi_MPV.append([low, high, MPV])
        bin_mid = 0.5 * (b["lo"] + b["hi"])
        if (not np.isnan(bin_mid)) and (not np.isnan(std68)) and n_jets_in_bin > n_jets_in_bin_min:
            bin_mid_points.append(bin_mid)
            sigmaEoverE.append(std68 / MPV if divide_by_MPV else std68)
            responses.append(MPV)
            y_ds, edges_ds = downsample_for_dashboard(y_normalized, edges)
            bin_records.append(
                {
                    "lo": b["lo"],
                    "hi": b["hi"],
                    "low": low,
                    "high": high,
                    "mpv": MPV,
                    "n_jets": n_jets_in_bin,
                    "edges": edges_ds,
                    "y": y_ds,
                }
            )
            print(
                f"Bin [{b['lo']}, {b['hi']}]: {sigma_method} = {std68:.4f}, low = {low:.4f}, "
                f"high = {high:.4f}, MPV={MPV}, N_in_bin={n_jets_in_bin}"
            )
        else:
            print(
                "NaN encountered in bin mid-point calculation / not enough statistics. "
                f"Number of jets in bin: {n_jets_in_bin}"
            )

    ax_hist[0].legend(fontsize=9)
    ax_hist[0].set_xlabel(x_label)
    ax_hist[0].set_ylabel("Entries")
    ax_hist[1].set_xlabel(x_label)
    ax_hist[1].set_ylabel("Entries")
    ax_hist[1].set_xlim([wmin, wmax])
    ax_hist[2].set_xlim([wmin, wmax])
    ax_hist[2].set_yscale("log")
    ax_hist[2].set_xlabel(x_label)
    ax_hist[2].set_ylabel("Entries")
    fig_hist.tight_layout()
    return (
        bin_mid_points,
        sigmaEoverE,
        fig_hist,
        responses,
        bins_to_histograms,
        lo_hi_MPV,
        field_names,
        total_statistics,
        bin_records,
    )


bin_to_histograms_storage = {}
bin_to_histograms_storage_neutral = {}
method_low_high_mid_point_storage = {}
fit_storage_Theta = {}
fit_storage_Phi = {}
fit_storage_Eta = {}
process_popt_storage = {}
dashboard_data = {p: {"energy": {}, "angles": {}} for p in processList}

method_to_color = {"std68": "blue", "gaussian_fit": "purple"}

jet_parts_to_process = ["_all"] if args.angles_only else ["_photons", "_neutral", "_charged", "_all"]
for jet_part in jet_parts_to_process:
    if not args.angles_only:
        fig_resolution_per_process, ax_resolution_per_process = plt.subplots(
            len(processList), 2, figsize=(8, 4 * len(processList)), sharex=False
        )
        fig_resolution_per_process_Njets, ax_resolution_per_process_Njets = plt.subplots(
            2, 2, figsize=(9, 9), sharex=False
        )
    for method in ["std68"]:
        print("-----------------------------------------------------------")
        print("Using peak width method:", method)
        if not args.angles_only:
            fig, ax = plt.subplots(
                2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]}
            )
        if jet_part == "_all":
            fig_theta, ax_theta = plt.subplots(
                2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]}
            )
            fig_phi, ax_phi = plt.subplots(
                2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]}
            )
            fig_eta, ax_eta = plt.subplots(
                2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]}
            )
        for proc_idx, process in enumerate(sorted(list(processList.keys()))):
            print("Process:", process)
            proc_data = raw["data"][process]
            if jet_part == "_all":
                (E_theta, sigma_theta, fig_theta_hist, response_theta, _, results_theta, _, total_stats_theta, bins_theta) = compute_resolution_for_process(
                    proc_data["angles"]["theta"],
                    sigma_method=method,
                    wmin=-0.05,
                    wmax=0.05,
                    divide_by_MPV=False,
                )
                xs_theta, ys_theta, popt_theta, pcov_theta = fit_resolution_model(
                    E_theta, sigma_theta, model="two_param"
                )
                fit_storage_Theta[process] = (
                    popt_theta, pcov_theta, xs_theta, ys_theta, E_theta, sigma_theta, results_theta,
                )
                if total_stats_theta >= 50000:
                    fig_theta_hist.savefig(os.path.join(outputDir, "bins_theta_{}.pdf".format(process)))
                else:
                    print(f"Skipping bins_theta plot for process {process}: N={total_stats_theta} < 50000")

                (E_phi, sigma_phi, fig_phi_hist, response_phi, _, results_phi, _, total_stats_phi, bins_phi) = compute_resolution_for_process(
                    proc_data["angles"]["phi"],
                    sigma_method=method,
                    wmin=-0.05,
                    wmax=0.05,
                    divide_by_MPV=False,
                    x_label="$\Delta \phi = \phi_{reco} - \phi_{true}$ [rad]",
                )
                if total_stats_phi >= 50000:
                    fig_phi_hist.savefig(os.path.join(outputDir, "bins_phi_{}.pdf".format(process)))
                else:
                    print(f"Skipping bins_phi plot for process {process}: N={total_stats_phi} < 50000")
                clr = PROCESS_COLORS.get(process, f"C{proc_idx}")
                xs_phi, ys_phi, popt_phi, pcov_phi = fit_resolution_model(
                    E_phi, sigma_phi, model="two_param"
                )
                fit_storage_Phi[process] = (
                    popt_phi, pcov_phi, xs_phi, ys_phi, E_phi, sigma_phi, results_phi,
                )
                ax_theta[0].plot(E_theta, sigma_theta, "x", color=clr)
                ax_theta[0].plot(
                    xs_theta, ys_theta, LINE_STYLES[process], color=clr,
                    label=HUMAN_READABLE_PROCESS_NAMES[process] + f" {print_params(popt_theta)}",
                )
                ax_theta[1].plot(E_theta, response_theta, ".--", label=process, color=clr)
                ax_phi[0].plot(E_phi, sigma_phi, "x", color=clr)
                ax_phi[0].plot(
                    xs_phi, ys_phi, LINE_STYLES[process], color=clr,
                    label=HUMAN_READABLE_PROCESS_NAMES[process] + f" {print_params(popt_phi)}",
                )
                ax_phi[1].plot(E_phi, response_phi, ".--", label=process, color=clr)
                ax_phi[1].set_xlabel("$E_{true}$ [GeV]")
                ax_phi[1].set_ylabel("Response in $\phi$")

                ax_theta[0].set_xlabel("$E_{true}$ [GeV]")
                ax_theta[0].set_ylabel(r"$\sigma_{\theta}$ [rad]")
                ax_theta[1].set_xlabel("$E_{true}$ [GeV]")
                ax_theta[1].set_ylabel("Response in $\\theta$")
                ax_phi[0].set_xlabel("$E_{true}$ [GeV]")
                ax_phi[0].set_ylabel(r"$\sigma_{\phi}$ [rad]")

                (E_eta, sigma_eta, fig_eta_hist, response_eta, _, results_eta, _, total_stats_eta, bins_eta_recs) = compute_resolution_for_process(
                    proc_data["angles"]["eta"],
                    sigma_method=method,
                    wmin=-0.05,
                    wmax=0.05,
                    divide_by_MPV=False,
                    x_label="$\Delta \eta = \eta_{reco} - \eta_{true}$",
                )
                if total_stats_eta >= 50000:
                    fig_eta_hist.savefig(os.path.join(outputDir, "bins_deltaEta_{}.pdf".format(process)))
                else:
                    print(f"Skipping bins_eta plot for process {process}: N={total_stats_eta} < 50000")
                xs_eta, ys_eta, popt_eta, pcov_eta = fit_resolution_model(
                    E_eta, sigma_eta, model="two_param"
                )
                fit_storage_Eta[process] = (
                    popt_eta, pcov_eta, xs_eta, ys_eta, E_eta, sigma_eta, results_eta,
                )
                ax_eta[0].plot(E_eta, sigma_eta, "x", color=clr)
                ax_eta[0].plot(
                    xs_eta, ys_eta, LINE_STYLES[process], color=clr,
                    label=HUMAN_READABLE_PROCESS_NAMES[process] + f" {print_params(popt_eta)}",
                )
                ax_eta[1].plot(E_eta, response_eta, ".--", label=process, color=clr)
                ax_eta[0].set_xlabel("$E_{true}$ [GeV]")
                ax_eta[0].set_ylabel(r"$\sigma_{\eta}$")
                ax_eta[1].set_xlabel("$E_{true}$ [GeV]")
                ax_eta[1].set_ylabel("Response in $\eta$")

                dashboard_data[process]["angles"]["theta"] = {
                    "mid_points": E_theta, "sigma": sigma_theta, "response": response_theta,
                    "model": "two_param", "popt": popt_theta.tolist(), "fit_x": xs_theta.tolist(), "fit_y": ys_theta.tolist(),
                    "bins": bins_theta,
                }
                dashboard_data[process]["angles"]["phi"] = {
                    "mid_points": E_phi, "sigma": sigma_phi, "response": response_phi,
                    "model": "two_param", "popt": popt_phi.tolist(), "fit_x": xs_phi.tolist(), "fit_y": ys_phi.tolist(),
                    "bins": bins_phi,
                }
                dashboard_data[process]["angles"]["eta"] = {
                    "mid_points": E_eta, "sigma": sigma_eta, "response": response_eta,
                    "model": "two_param", "popt": popt_eta.tolist(), "fit_x": xs_eta.tolist(), "fit_y": ys_eta.tolist(),
                    "bins": bins_eta_recs,
                }
            if args.angles_only:
                continue
            (
                bin_mid_points, sigmaEoverE, fig_histograms, resp, bin_to_histograms,
                mpv_lo_hi, field_names, _, bin_records,
            ) = compute_resolution_for_process(
                proc_data["energy"][jet_part],
                sigma_method=method,
            )
            if process not in method_low_high_mid_point_storage:
                method_low_high_mid_point_storage[process] = {}
            if jet_part == "_all":
                method_low_high_mid_point_storage[process][method] = mpv_lo_hi
            if method == "std68" and jet_part == "_all":
                bin_to_histograms_storage[process] = bin_to_histograms
                fig_histograms.tight_layout()
                fig_histograms.savefig(os.path.join(outputDir, "bins_{}_{}.pdf".format(process, method)))
            if method == "std68" and jet_part == "_neutral":
                bin_to_histograms_storage_neutral[process] = bin_to_histograms
                fig_histograms.tight_layout()
                fig_histograms.savefig(os.path.join(outputDir, "bins_NEUTRAL_{}_{}.pdf".format(process, method)))
            clr = PROCESS_COLORS.get(process, f"C{proc_idx}")
            if len(bin_mid_points) < 2:
                print(f"Not enough points to fit for process {process} using method {method}. Skipping.")
                continue
            xs, ys, popt, pcov = fit_resolution_model(
                bin_mid_points, sigmaEoverE, model="three_param",
                bounds_override=([0.0, 0.005, 0.0], [np.inf, 0.04, np.inf]),
            )
            print(f"Fitted parameters for {process} using {method}: {popt}")
            if process not in process_popt_storage:
                process_popt_storage[process] = {}
            process_popt_storage[process][method + jet_part] = (
                popt, pcov, xs, ys, bin_mid_points, sigmaEoverE, mpv_lo_hi, field_names,
            )
            dashboard_data[process]["energy"][jet_part] = {
                "mid_points": bin_mid_points, "sigma_over_E": sigmaEoverE, "response": resp,
                "model": "three_param", "popt": popt.tolist(), "fit_x": xs.tolist(), "fit_y": ys.tolist(),
                "bins": bin_records,
            }
            ax[0].plot(bin_mid_points, sigmaEoverE, "x", color=clr)
            ax[0].plot(
                xs, ys, LINE_STYLES[process], color=clr,
                label=HUMAN_READABLE_PROCESS_NAMES[process] + f" {print_params(popt)}",
            )
            ax[1].plot(bin_mid_points, resp, ".--", label=process, color=clr)
            if method in method_to_color:
                ax_resolution_per_process[proc_idx, 0].plot(
                    bin_mid_points, sigmaEoverE, "x", label=method + f" {print_params(popt)}",
                    color=method_to_color[method],
                )
                ax_resolution_per_process[proc_idx, 0].plot(xs, ys, LINE_STYLES[process], color=method_to_color[method])
                ax_resolution_per_process[proc_idx, 0].set_title(HUMAN_READABLE_PROCESS_NAMES[process])
                ax_resolution_per_process[proc_idx, 1].plot(
                    bin_mid_points, resp, ".--", label=method, color=method_to_color[method],
                )
            if method == "std68":
                if LINE_STYLES.get(process, "") == "-":
                    row = 0
                elif LINE_STYLES.get(process, "") == ":":
                    row = 1
                else:
                    row = None
                if row is not None:
                    ax_resolution_per_process_Njets[row, 0].plot(
                        bin_mid_points, sigmaEoverE, "x",
                        label=HUMAN_READABLE_PROCESS_NAMES[process] + f" {print_params(popt)}", color=clr,
                    )
                    ax_resolution_per_process_Njets[row, 0].plot(xs, ys, "--", color=clr)
                    ax_resolution_per_process_Njets[row, 1].plot(
                        bin_mid_points, resp, ".--", label=HUMAN_READABLE_PROCESS_NAMES[process], color=clr,
                    )
            if not args.angles_only:
                ax_resolution_per_process[proc_idx, 0].set_xlabel("$E_{true}$ [GeV]")
                ax_resolution_per_process[proc_idx, 0].set_ylabel(r"$\sigma_E / E$")
                ax_resolution_per_process[proc_idx, 1].set_xlabel("$E_{true}$ [GeV]")
                ax_resolution_per_process[proc_idx, 1].set_ylabel("Response")
                ax_resolution_per_process[proc_idx, 0].legend()
                ax_resolution_per_process[proc_idx, 0].grid(True)
                ax_resolution_per_process[proc_idx, 1].grid(True)
        if not args.angles_only:
            ax_resolution_per_process_Njets[0, 0].set_title("Final state containing b-jets")
            ax_resolution_per_process_Njets[1, 0].set_title("Final state containing only light and gluon jets")
            ax_resolution_per_process_Njets[0, 0].legend(title="l ∈ {u, d, s}; q ∈ {u, d, s, c, b}", fontsize=9.5, title_fontsize=8)
            ax_resolution_per_process_Njets[1, 0].legend(title="l ∈ {u, d, s}; q ∈ {u, d, s, c, b}", fontsize=9.5, title_fontsize=8)
            ax_resolution_per_process_Njets[0, 0].set_xlabel("$E_{true}$ [GeV]")
            ax_resolution_per_process_Njets[0, 1].set_xlabel("$E_{true}$ [GeV]")
            ax_resolution_per_process_Njets[1, 0].set_xlabel("$E_{true}$ [GeV]")
            ax_resolution_per_process_Njets[1, 1].set_xlabel("$E_{true}$ [GeV]")
            ax_resolution_per_process_Njets[0, 0].set_ylabel(r"$\sigma_E / E$")
            ax_resolution_per_process_Njets[1, 0].set_ylabel(r"$\sigma_E / E$")
            ax_resolution_per_process_Njets[1, 1].set_ylabel(r"Response in E reco / E true")
            ax_resolution_per_process_Njets[0, 1].set_ylabel(r"Response in E reco / E true")
            ax_resolution_per_process_Njets[0, 0].grid()
            ax_resolution_per_process_Njets[1, 1].grid()
            ax_resolution_per_process_Njets[0, 1].grid()
            ax_resolution_per_process_Njets[1, 0].grid()
            ax[0].legend()
            ax[0].set_xlabel("$E_{true}$ [GeV]")
            ax[0].set_ylabel(r"$\sigma_E / E$")
            ax[0].set_title(r"Jet Energy Resolution ($\frac{A}{\sqrt{E}}$ ⊕ \frac{B}{E} ⊕ $C$)")
            ax[0].grid(True, alpha=0.3)
            ax[1].set_ylabel("Response")
            ax[1].set_xlabel("$E_{true}$ [GeV]")
            ax[1].grid()
            fig.tight_layout()
            fig.savefig(os.path.join(outputDir, "jet_energy_resolution_{}.pdf".format(method)))
        if jet_part == "_all":
            ax_theta[0].legend()
            ax_theta[0].grid(True)
            ax_theta[1].grid(True)
            fig_theta.tight_layout()
            fig_theta.savefig(os.path.join(outputDir, "jet_theta_resolution_{}.pdf".format(method)))
            ax_phi[0].legend()
            ax_phi[0].grid(True)
            ax_phi[1].grid(True)
            fig_phi.tight_layout()
            fig_phi.savefig(os.path.join(outputDir, "jet_phi_resolution_{}.pdf".format(method)))
            ax_eta[0].legend()
            ax_eta[0].grid(True)
            ax_eta[1].grid(True)
            fig_eta.tight_layout()
            fig_eta.savefig(os.path.join(outputDir, "jet_eta_resolution_{}.pdf".format(method)))
    if not args.angles_only:
        fig_resolution_per_process.tight_layout()
        fig_resolution_per_process_Njets.tight_layout()
        fig_resolution_per_process.savefig(
            os.path.join(outputDir, "jet_energy_resolution_per_process_comparison{}.pdf".format(jet_part))
        )
        fig_resolution_per_process_Njets.savefig(
            os.path.join(outputDir, "jet_energy_resolution_per_process_comparison_Njets{}.pdf".format(jet_part))
        )

pickle.dump(
    {"theta": fit_storage_Theta, "phi": fit_storage_Phi, "eta": fit_storage_Eta},
    open(os.path.join(outputDir, "angle_fit_params_per_process.pkl"), "wb"),
)

if not args.angles_only:
    pickle.dump(
        process_popt_storage,
        open(os.path.join(outputDir, "energy_fit_params_per_process.pkl"), "wb"),
    )
    method_to_color = {
        "std68": "blue", "RMS": "orange", "interquantile_range": "green", "DSCB": "red", "gaussian_fit": "purple",
    }
else:
    pickle.dump(dashboard_data, open(os.path.join(outputDir, "resolution_dashboard_data.pkl"), "wb"))
    print("Saved dashboard data to:", os.path.join(outputDir, "resolution_dashboard_data.pkl"))
    import sys
    sys.exit(0)

### Plot each bin on a separate plot, but different processes on same plot ###
fig, ax = plt.subplots(len(binsE) - 1, 1, figsize=(6, 4 * (len(binsE) - 1)), sharex=False)
fig_bins, ax_bins = plt.subplots(len(binsE) - 1, 1, figsize=(6, 4 * (len(binsE) - 1)), sharex=False)

for i in range(len(binsE) - 1):
    for process in sorted(list(processList.keys())):
        y_normalized, edges = bin_to_histograms_storage[process][i]
        ax[i].step(edges[:-1], y_normalized, where="post", label=process)
        ax_bins[i].step(edges[:-1], y_normalized, where="post", label=process)
        for method in method_low_high_mid_point_storage[process]:
            lo, hi, mpv = method_low_high_mid_point_storage[process][method][i]
            for _ax in [ax[i], ax_bins[i]]:
                _ax.axvline(lo, color=method_to_color[method], linestyle="--", alpha=0.8)
                _ax.axvline(hi, color=method_to_color[method], linestyle="--", alpha=0.8)
                _ax.axvline(mpv, color=method_to_color[method], linestyle="-", alpha=0.8)
    ax[i].set_title(f"Bin [{binsE[i]}, {binsE[i + 1]}] GeV")
    ax[i].set_ylabel("Entries")
    ax[i].legend()
    ax[i].set_xlim([0.95, 1.05])
    ax_bins[i].set_title(f"Bin [{binsE[i]}, {binsE[i + 1]}] GeV")
    ax_bins[i].set_ylabel("Entries")
    ax_bins[i].legend()
    ax_bins[i].set_yscale("log")

ax[-1].set_xlabel(r"$E_{reco} / E_{true}$")
ax_bins[-1].set_xlabel(r"$E_{reco} / E_{true}$")

fig.tight_layout()
fig.savefig(os.path.join(outputDir, "jet_energy_bins_comparison.pdf"))

fig_bins.tight_layout()
fig_bins.savefig(os.path.join(outputDir, "jet_energy_bins_comparison_full_axis.pdf"))

eta_scan_data = {}
for method in ["std68"]:
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]})
    for process in sorted(list(processList.keys())):
        bin_mid_points, sigmaEoverE, fig_histograms, resp, _, _, _, _, bin_records = compute_resolution_for_process(
            raw["data"][process]["eta_binned_E"], sigma_method=method,
        )
        if method == "std68":
            fig_histograms.tight_layout()
            fig_histograms.savefig(os.path.join(outputDir, "bins_E_vs_eta_{}.pdf".format(process)))
        ax[0].plot(bin_mid_points, sigmaEoverE, LINE_STYLES[process], label=HUMAN_READABLE_PROCESS_NAMES[process], color=PROCESS_COLORS[process])
        ax[1].plot(bin_mid_points, resp, LINE_STYLES[process], label=HUMAN_READABLE_PROCESS_NAMES[process], color=PROCESS_COLORS[process])
        eta_scan_data[process] = {"mid_points": bin_mid_points, "sigma_over_E": sigmaEoverE, "response": resp, "bins": bin_records}
    ax[0].legend()
    ax[0].set_xlabel("$\eta$")
    ax[0].set_ylabel(r"$\sigma_E / E$")
    ax[0].set_title("Jet Energy Resolution vs $\Theta$")
    ax[0].grid(True, alpha=0.3)
    ax[1].grid(True, alpha=0.3)
    ax[1].set_title("Energy Response vs $\Theta$")
    ax[1].set_xlabel("$\eta$")
    ax[1].set_ylabel("Response")
    ax[1].grid()
    fig.tight_layout()
    fig.savefig(os.path.join(outputDir, "JER_vs_eta_{}.pdf".format(method)))

costheta_scan_data = {}
for method in ["std68"]:
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]})
    for process in sorted(list(processList.keys())):
        bin_mid_points, sigmaEoverE, fig_histograms, resp, _, _, _, _, bin_records = compute_resolution_for_process(
            raw["data"][process]["costheta_binned_E"], sigma_method=method,
        )
        if method == "std68":
            fig_histograms.tight_layout()
            fig_histograms.savefig(os.path.join(outputDir, "bins_E_vs_CosTheta_{}.pdf".format(process)))
        ax[0].plot(bin_mid_points, sigmaEoverE, "x", color=PROCESS_COLORS[process])
        ax[1].plot(bin_mid_points, resp, "x", color=PROCESS_COLORS[process])
        ax[0].plot(bin_mid_points, sigmaEoverE, LINE_STYLES[process], label=HUMAN_READABLE_PROCESS_NAMES[process], color=PROCESS_COLORS[process])
        ax[1].plot(bin_mid_points, resp, LINE_STYLES[process], label=HUMAN_READABLE_PROCESS_NAMES[process], color=PROCESS_COLORS[process])
        costheta_scan_data[process] = {"mid_points": bin_mid_points, "sigma_over_E": sigmaEoverE, "response": resp, "bins": bin_records}
    ax[0].legend()
    ax[0].set_xlabel(r"cos $\theta$ [GeV]")
    ax[0].set_ylabel(r"$\sigma_E / E$")
    ax[0].set_title("Jet Energy Resolution vs Jet Energy")
    ax[0].grid(True, alpha=0.3)
    ax[1].grid(True, alpha=0.3)
    ax[1].set_title("Energy Resolution vs Jet Angle")
    ax[1].set_xlabel(r"cos $\theta$ [GeV]")
    ax[1].set_ylabel("Response")
    ax[1].grid()
    fig.tight_layout()
    fig.savefig(os.path.join(outputDir, "JER_vs_CosTheta_resolution_{}.pdf".format(method)))

for process in dashboard_data:
    dashboard_data[process]["eta_scan"] = eta_scan_data.get(process)
    dashboard_data[process]["costheta_scan"] = costheta_scan_data.get(process)

pickle.dump(dashboard_data, open(os.path.join(outputDir, "resolution_dashboard_data.pkl"), "wb"))
print("Saved dashboard data to:", os.path.join(outputDir, "resolution_dashboard_data.pkl"))
