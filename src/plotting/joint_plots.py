import matplotlib.pyplot as plt
import numpy as np
import pickle
from ROOT import TFile
from scipy.optimize import curve_fit
import argparse
import os
import sys
from src.process_config import (
    PROCESS_COLORS,
    HUMAN_READABLE_PROCESS_NAMES,
    LINE_STYLES,
    NUMBER_OF_JETS,
    PROCESS_TO_ROW_COL,
)

parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", type=str, required=True)
parser.add_argument(
    "--AK-comparison", action="store_true"
)  # If turned on, it will produce plots for comparison of anti-kt and Durham
parser.add_argument(
    "--energy-recovery", action="store_true"  # If toggled, it will turn on the
)
args = parser.parse_args()

inputDir = args.inputDir
outputDir = os.path.join(inputDir, "plots")
if args.AK_comparison:
    if args.energy_recovery:
        outputDir = os.path.join(outputDir, "comparison_AK_energy_recovery")
    else:
        outputDir = os.path.join(outputDir, "comparison_AK")
os.makedirs(outputDir, exist_ok=True)

bins_E = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]


def annotate_matrix_plot_with_arrows(fig):
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    ax.set_axis_off()
    ax.annotate(
        "",
        xy=(0.009, 0.60),
        xycoords="figure fraction",
        xytext=(0.009, 0.99),
        textcoords="figure fraction",
        arrowprops=dict(arrowstyle="->", lw=1.2),
        color="gray",
    )
    ax.text(
        0.0048,
        0.75,
        "More B-hadron content",
        transform=ax.transAxes,
        rotation=90,
        ha="center",
        va="center",
        fontsize=9.5,
    )
    ax.annotate(
        "",
        xy=(0.40, 0.99),
        xycoords="figure fraction",
        xytext=(0.01, 0.99),
        textcoords="figure fraction",
        arrowprops=dict(arrowstyle="->", lw=1.2),
        color="gray",
    )
    ax.text(
        0.25,
        0.995,
        "Higher number of final-state jets",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=9.5,
    )


def print_params(popt):
    if len(popt) == 2:
        return f"S={round(popt[0], 2)} C={round(popt[1], 2)}"
    return f"S={round(popt[0], 2)} N={round(popt[2], 2)} C={round(popt[1], 2)}"


def print_noise_model(popt):
    # print the noise model based on the number of parameters (2: A/sqrt(E) + C ; otherwise, A/sqrt(E) + B + C/E)
    if len(popt) == 2:
        return r"$\frac{A}{\sqrt{E}} \oplus C$"
    return r"$\frac{A}{\sqrt{E}} \oplus B \oplus \frac{C}{E}$"


def root_file_get_hist_and_edges(root_file, hist_name, rebin_factor=1):
    # Assume root file is a text, open it
    root_file_obj = TFile.Open(root_file, "READ")
    h = root_file_obj.Get(hist_name)
    if not h:
        raise Exception(f"Warning: histogram {hist_name} not found")
    nb = h.GetNbinsX()
    edges = np.array(
        [h.GetXaxis().GetBinLowEdge(1)]
        + [h.GetXaxis().GetBinUpEdge(i) for i in range(1, nb + 1)]
    )
    y = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    assert len(edges) == len(y) + 1
    return y, edges


if args.energy_recovery:
    AK_prefix = "PF_E_recovery_AntiKtR"
else:
    AK_prefix = "PF_AntiKtR"

if args.AK_comparison:
    method_dict = {"PF_Durham": "Durham"}
    for radius in [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]:
        radius_str = int(radius * 10)
        if len(str(radius_str)) == 1:
            radius_str = f"0{radius_str}"
        if args.energy_recovery:
            method_dict[f"{AK_prefix}{radius_str}"] = f"AK{radius}-ER"
        else:
            method_dict[f"{AK_prefix}{radius_str}"] = f"AK{radius}"
else:
    method_dict = {
        "PF_Durham_IdealMatching": "PFJets + Ideal Matching",
        "PF_Durham": "PFJets",
        "CaloJets_Durham": "CaloJets",
    }

process_for_detailed_bins_plots = [
    "p8_ee_ZH_vvbb_ecm240",
    "p8_ee_ZH_vvgg_ecm240",
    "p8_ee_ZH_vvqq_ecm240",
]

twojet_processes = [
    "p8_ee_ZH_vvbb_ecm240",
    "p8_ee_ZH_vvgg_ecm240",
    "p8_ee_ZH_vvqq_ecm240",
]

gluon_comparison_processes = {  # process idx to plot idx - group together by 2 jets and 4 jets
    # Two jets
    "p8_ee_ZH_vvgg_ecm240": 0,
    "p8_ee_ZH_vvqq_ecm240": 0,
    # Four jets
    "p8_ee_ZH_qqqq_ecm240": 1,
    "p8_ee_ZH_qqgg_ecm240": 1,
    "p8_ee_ZH_bbgg_ecm240": 1,
    "p8_ee_ZH_bbbb_ecm240": 1,
}

# Resolution plots, PF vs PF + Ideal matching
fig, ax = plt.subplots(5, 3, figsize=(10, 11))
fig_mH, ax_mH = plt.subplots(2, 3, figsize=(9, 7))
fig_mH_per_process, ax_mH_per_process = plt.subplots(5, 3, figsize=(10, 11))
fig_mH_twojets, ax_mH_twojets = plt.subplots(4, 2, figsize=(9, 9))
fig_E_mH_gluons, ax_E_mH_gluons = plt.subplots(2, 3, figsize=(12, 8))


def get_func_fit(
    mid_points,
    Rs,
    confusion_term=True,
    min_E=0.0,
    bounds_constant=[0, np.inf],
    bounds_confusion=[0, np.inf],
    bounds_noise=[0, np.inf],
    error_summation=False,
):
    # basically filter out points with mid_points < min_E
    mask = np.array(mid_points) >= min_E
    mid_points = np.array(mid_points)[mask]
    Rs = np.array(Rs)[mask]
    if confusion_term:

        def resolution_func(E, a, b, c):
            if error_summation:
                return np.sqrt((a / np.sqrt(E)) ** 2 + b**2 + (c / E) ** 2)
            return a / np.sqrt(E) + b + c / E

        bounds = (
            [bounds_noise[0], bounds_constant[0], bounds_confusion[0]],
            [bounds_noise[1], bounds_constant[1], bounds_confusion[1]],
        )
        popt, pcov = curve_fit(
            resolution_func,
            mid_points,
            Rs,
            p0=[0.5, 0.03, 0.1],
            maxfev=10000,
            bounds=bounds,
        )
        # return a smooth function throughout mid_points
    else:

        def resolution_func(E, a, b):
            if error_summation:
                return np.sqrt((a / np.sqrt(E)) ** 2 + b**2)
            return a / np.sqrt(E) + b

        bounds = (
            [bounds_noise[0], bounds_constant[0]],
            [bounds_noise[1], bounds_constant[1]],
        )
        popt, pcov = curve_fit(
            resolution_func, mid_points, Rs, p0=[0.5, 0.03], maxfev=10000, bounds=bounds
        )
        # return a smooth function throughout mid_points
    xs = np.linspace(min(mid_points), max(mid_points), 100)
    ys = resolution_func(xs, *popt)
    popt = np.abs(popt)
    return xs, ys, popt, pcov


if args.AK_comparison:
    methods_filtered = sorted(method_dict.keys())
    method_color = {"PF_Durham": "#1f77b4"}  # nice blue for Durham
    # build purple shades for AntiKt radii, getting darker with radius
    ak_keys = [k for k in methods_filtered if k.startswith(AK_prefix)]

    def _radius_key(k):
        try:
            return int(k.split("R")[-1])
        except Exception:
            return 0

    ak_keys_sorted = sorted(ak_keys, key=_radius_key)
    cmap = plt.cm.Purples
    vals = np.linspace(0.4, 0.9, max(1, len(ak_keys_sorted)))

    def _rgba_to_hex(rgba):
        r, g, b, _ = rgba
        return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))

    for i, k in enumerate(ak_keys_sorted):
        method_color[k] = _rgba_to_hex(cmap(vals[i]))
    method_linestyle = {key: "-" for key in methods_filtered}

else:
    methods_filtered = ["PF_Durham", "PF_Durham_IdealMatching"]
    method_color = {
        "PF_Durham": "blue",
        "PF_Durham_IdealMatching": "orange",
    }
    method_linestyle = {  # For the m_H plots
        "PF_Durham": "-",
        "PF_Durham_IdealMatching": "--",
    }
for method in methods_filtered:
    f = pickle.load(
        open(
            os.path.join(
                inputDir, f"{method}/plots_resolution/energy_fit_params_per_process.pkl"
            ),
            "rb",
        )
    )
    f_mH = pickle.load(
        open(
            os.path.join(
                inputDir, f"{method}/plots_mass/Higgs_mass_histograms_data.pkl"
            ),
            "rb",
        )
    )
    for i, process in enumerate(PROCESS_COLORS.keys()):
        if process not in f:
            continue
        color = PROCESS_COLORS[process]
        label = HUMAN_READABLE_PROCESS_NAMES[process]
        linestyle = "--"
        cp = f[process]["std68_all"]
        xs, ys = cp[2], cp[3]
        x_pts, y_pts = cp[4], cp[5]
        col, row = PROCESS_TO_ROW_COL[process]
        # get func fit
        xs, ys, fit_params, fit_cov = get_func_fit(x_pts, y_pts, confusion_term=False)
        ax[row, col].plot(
            xs,
            ys,
            label=f"{method_dict[method]} {print_params(fit_params)}",
            color=method_color[method],
            linestyle=linestyle,
        )
        ax[row, col].plot(x_pts, y_pts, "x", color=method_color[method], markersize=4)
        ax[row, col].set_ylabel("$\sigma_E / E_{true}$")
        # Set E_true as xlabel
        ax[row, col].set_xlabel("$E_{true}$ [GeV]")
        ax[row, col].set_title(label)  # + f" ({print_noise_model(fit_params)})")
        Higgs_x, Higgs_y = f_mH[process]["x_vals_reco"], f_mH[process]["y_vals_reco"]
        Higgs_bin_width = Higgs_x[1] - Higgs_x[0] if len(Higgs_x) > 1 else 1.0
        Higgs_edges = np.concatenate(
            (Higgs_x - Higgs_bin_width / 2, [Higgs_x[-1] + Higgs_bin_width / 2])
        )
        if process in gluon_comparison_processes and method == "PF_Durham":
            plot_idx = gluon_comparison_processes[process]
            ax_E_mH_gluons[plot_idx, 1].hist(
                Higgs_x,
                bins=Higgs_edges,
                weights=Higgs_y,
                histtype="step",
                label=f"{label} ({method_dict[method]})",
            )
            ax_E_mH_gluons[plot_idx, 2].hist(
                Higgs_x,
                bins=Higgs_edges,
                weights=Higgs_y,
                histtype="step",
                label=f"{label} ({method_dict[method]})",
            )
            # Set axes for the plots 1 and 2
            ax_E_mH_gluons[plot_idx, 2].set_xlabel("$m_H$ (reco.) [GeV]")
            ax_E_mH_gluons[plot_idx, 1].set_xlabel("$m_H$ (reco.) [GeV]")
            ax_E_mH_gluons[plot_idx, 1].set_xlabel("Normalized events")
            ax_E_mH_gluons[plot_idx, 2].set_xlabel("Normalized events")
            n_jets = [2, 4][plot_idx]
            ax_E_mH_gluons[plot_idx, 1].set_title(f"{n_jets} jets; $m_H$")
            ax_E_mH_gluons[plot_idx, 2].set_title(f"{n_jets} jets; $m_H$ (zoom)")
            ax_E_mH_gluons[plot_idx, 0].set_title(
                f"{n_jets} jets; Jet Energy Resolution"
            )
            ax_E_mH_gluons[plot_idx, 2].set_xlim(115, 135)
            clr = "C{}".format(i)
            ax_E_mH_gluons[plot_idx, 0].plot(
                xs,
                ys,
                label=f"{label} {print_params(fit_params)}",
                linestyle=linestyle,
                color=clr,
            )
            ax_E_mH_gluons[plot_idx, 0].plot(x_pts, y_pts, "x", markersize=4, color=clr)
            ax_E_mH_gluons[plot_idx, 0].set_ylabel("$\sigma_E / E_{true}$")
            ax_E_mH_gluons[plot_idx, 0].set_xlabel("$E_{true}$ [GeV]")
            # Also plot the energy resolution histogram on (plot_idx, 0)
        if LINE_STYLES[process] == ":":
            for a in ax_mH[0]:
                a.hist(
                    Higgs_x,
                    bins=Higgs_edges,
                    weights=Higgs_y,
                    histtype="step",
                    label=f"{label} ({method_dict[method]})",
                    color=PROCESS_COLORS[process],
                    linestyle=method_linestyle[method],
                )
            ax_mH[0, 0].set_title(r"H → Light jets")
        elif LINE_STYLES[process] == "-":
            for a in ax_mH[1]:
                a.hist(
                    Higgs_x,
                    bins=Higgs_edges,
                    weights=Higgs_y,
                    histtype="step",
                    label=f"{label} ({method_dict[method]})",
                    color=PROCESS_COLORS[process],
                    linestyle=method_linestyle[method],
                )
            ax_mH[1, 0].set_title(r"H → Light and b-jets")
        ax_mH_per_process[row, col].hist(
            Higgs_x,
            bins=Higgs_edges,
            weights=Higgs_y,
            histtype="step",
            label=f"{method_dict[method]}",
            color=method_color[method],
            linestyle=method_linestyle[method],
        )
        ax_mH_per_process[row, col].set_title(label)
        if NUMBER_OF_JETS.get(process) == 2 and method == "PF_Durham":
            ax_mH_twojets[0, 0].hist(
                Higgs_x,
                bins=Higgs_edges,
                weights=Higgs_y,
                histtype="step",
                label=f"{label} ({method_dict[method]})",
                linestyle=method_linestyle[method],
            )
            ax_mH_twojets[0, 1].hist(
                Higgs_x,
                bins=Higgs_edges,
                weights=Higgs_y,
                histtype="step",
                label=f"{label} ({method_dict[method]})",
                linestyle=method_linestyle[method],
            )
            ##########################
            ax_mH_twojets[0, 0].set_title(r"H → 2 jets")
            ax_mH_twojets[0, 1].set_title(r"H → 2 jets (zoom)")
            # If process is in twojet_processes, also plot it in the separate plot
            proc_idx = twojet_processes.index(process)
            # plot on the proc_idx +1  - th row in the ax_mH_twojets
            # get h_mH_reco from ../../idea_fullsim/fast_sim/{histograms_folder}/PFDurham_ISR/{process}.root
            root_file = os.path.join(inputDir, f"PF_Durham/{process}.root")
            y_reco_f, edges_reco_f = root_file_get_hist_and_edges(
                root_file, "h_mH_reco"
            )
            bin_centers_gen = 0.5 * (edges_reco_f[:-1] + edges_reco_f[1:])
            bin_centers_reco = 0.5 * (edges_reco_f[:-1] + edges_reco_f[1:])
            ax_mH_twojets[proc_idx + 1, 1].step(
                bin_centers_reco,
                y_reco_f,
                where="mid",
                label="PFJets",
                linestyle="solid",
                color="red",
                alpha=0.7,
            )
            ax_mH_twojets[proc_idx + 1, 0].step(
                bin_centers_reco,
                y_reco_f,
                where="mid",
                label="PFJets",
                linestyle="solid",
                color="red",
                alpha=0.7,
            )
            ax_mH_twojets[proc_idx + 1, 0].set_title(f"{label}")
            ax_mH_twojets[proc_idx + 1, 1].set_title(f"{label}")
            root_file = os.path.join(inputDir, f"PFDurham_ISR_NoFilter/{process}.root")
            if os.path.exists(root_file):
                y_reco, edges_reco = root_file_get_hist_and_edges(
                    root_file, "hist_calo_hist_E"
                )
                bin_centers_gen = 0.5 * (edges_reco[:-1] + edges_reco[1:])
                bin_centers_reco = 0.5 * (edges_reco[:-1] + edges_reco[1:])
                ax_mH_twojets[proc_idx + 1, 1].step(
                    bin_centers_reco,
                    y_reco,
                    where="mid",
                    label="Using all particles",
                    linestyle="solid",
                    color="black",
                    alpha=0.7,
                )
                ax_mH_twojets[proc_idx + 1, 0].step(
                    bin_centers_reco,
                    y_reco,
                    where="mid",
                    label="Using all particles",
                    linestyle="solid",
                    color="black",
                    alpha=0.7,
                )

for i in range(len(ax)):
    for j in range(len(ax[i])):
        ax[i, j].legend(fontsize=6.5)
        ax[i, j].grid()

for i in range(len(ax_E_mH_gluons)):
    for j in range(len(ax_E_mH_gluons[i])):
        ax_E_mH_gluons[i, j].legend(fontsize=7)
        ax_E_mH_gluons[i, j].grid()

for i in range(len(ax_mH)):
    for j in range(len(ax_mH[i])):
        ax_mH[i, j].set_xlabel("$m_H$ (reco.) [GeV]")
        ax_mH[i, j].set_xlabel("$m_H$ (reco.) [GeV]")
        ax_mH[i, j].grid()
    ax_mH[i, 2].set_xlim(60, 180)
    ax_mH[i, 2].set_yscale("log")
    ax_mH[i, 1].set_xlim(100, 140)
    ax_mH[i, 0].legend(title="q ∈ {u, d, s}", fontsize=7.5, title_fontsize=8)

for j in range(len(ax_mH_twojets)):
    for i in range(len(ax_mH_twojets[0])):
        ax_mH_twojets[j, i].set_xlabel("$m_H$ (reco.) [GeV]")
        ax_mH_twojets[j, i].grid()
        ax_mH_twojets[j, i].set_xlim(100, 140)
        ax_mH_twojets[j, i].legend(
            title="2-jet processes", fontsize=7.5, title_fontsize=8
        )

for j in range(len(ax_mH_per_process)):
    for i in range(len(ax_mH_per_process[j])):
        ax_mH_per_process[j, i].set_xlabel("$m_H$ (reco.) [GeV]")
        ax_mH_per_process[j, i].grid()
        ax_mH_per_process[j, i].set_xlim(90, 150)
        ax_mH_per_process[j, i].legend(fontsize=6.5)
fig.tight_layout()
annotate_matrix_plot_with_arrows(fig)
fig.tight_layout()
fig_mH_path_per_process = os.path.join(outputDir, f"Higgs_mass_per_process.pdf")

if args.AK_comparison:
    fig_path = os.path.join(outputDir, f"Jet_Energy_Resolution.pdf")
    fig_mH_path = os.path.join(outputDir, f"Higgs_mass.pdf")
    fig_mH_twojets_path = os.path.join(outputDir, f"Higgs_mass_2jets.pdf")
else:
    fig_path = os.path.join(outputDir, f"JER_comparison_PF_and_ideal_matching.pdf")
    fig_mH_path = os.path.join(
        outputDir, f"Higgs_mass_comparison_PF_and_ideal_matching.pdf"
    )
    fig_mH_twojets_path = os.path.join(
        outputDir, f"Higgs_mass_comparison_PF_and_ideal_matching_2jets.pdf"
    )

fig_mH.tight_layout()
fig_mH_twojets.tight_layout()
fig_mH_per_process.tight_layout()

print("Saving figure to", fig_mH_path)
print("Saving figure to", fig_mH_twojets_path)
print("Saving figure to", fig_path)
print("Saving figure to", fig_mH_path_per_process)

fig.savefig(fig_path)
fig_mH.savefig(fig_mH_path)
fig_mH_twojets.savefig(fig_mH_twojets_path)
fig_mH_per_process.savefig(fig_mH_path_per_process)

if args.AK_comparison:
    sys.exit(0)

#############################################
# Resolution plots, PF vs PF+Ideal matching #
#############################################

fig, ax = plt.subplots(5, 3, figsize=(10, 11))
fig_mH, ax_mH = plt.subplots(2, 2, figsize=(8, 6))
################################################################################################################
fig_fit_trials, ax_fit_trials = plt.subplots(5, 3, figsize=(10, 11))
fig_fit_trials_calojets, ax_fit_trials_calojets = plt.subplots(5, 3, figsize=(10, 11))
################################################################################################################

figs_processes, axs_processes = {}, {}
histogram_limits = [[0, 2], [0.3, 1.4], [0.75, 1.15], [0.95, 1.05], [0.99, 1.01]]
for method in [
    "PF_Durham",
    "CaloJets_Durham",
    "PF_Durham_IdealMatching",
]:
    method_color = {
        "PF_Durham": "blue",
        "CaloJets_Durham": "green",
        "PF_Durham_IdealMatching": "orange",
    }
    method_linestyle = {  # For the m_H plots
        "PF_Durham": "-",
        "CaloJets_Durham": "--",
        "PF_Durham_IdealMatching": ":",
    }
    f = pickle.load(
        open(
            os.path.join(
                inputDir,
                f"{method}/plots_resolution/energy_fit_params_per_process.pkl",
            ),
            "rb",
        )
    )
    f_mH = pickle.load(
        open(
            os.path.join(
                inputDir, f"{method}/plots_mass/Higgs_mass_histograms_data.pkl"
            ),
            "rb",
        )
    )
    for i, process in enumerate(PROCESS_COLORS.keys()):
        if process not in f:
            continue
        color = PROCESS_COLORS[process]
        label = HUMAN_READABLE_PROCESS_NAMES[process]
        if process in process_for_detailed_bins_plots:
            # Make detailed binning plots for these processes
            if method == "CaloJets_Durham":
                binning_metadata = f[process]["std68_all"]
            else:
                binning_metadata = f[process]["std68_all"]
            hist_names = binning_metadata[7]
            if process not in figs_processes:
                n_bins = len(binning_metadata[6])
                figs_processes[process], axs_processes[process] = plt.subplots(
                    n_bins,
                    len(histogram_limits),
                    figsize=(15, len(histogram_limits) * n_bins),
                )
            for b in range(n_bins):
                low_point, high_point, mpv_point = binning_metadata[6][b]
                y_hist, edges = root_file_get_hist_and_edges(
                    os.path.join(inputDir, f"{method}/{process}.root"),
                    hist_names[b],
                )
                # Now plot this in each column of bin b.
                # Each plot has different x axis limits (historam_limits)
                for col in range(len(histogram_limits)):
                    mask = (edges[:-1] >= histogram_limits[col][0]) & (
                        edges[1:] <= histogram_limits[col][1]
                    )
                    # Do a step histogram using mask
                    axs_processes[process][b, col].step(
                        edges[:-1][mask],
                        y_hist[mask],
                        where="post",
                        color=method_color[method],
                        label=f"{method_dict[method]}",
                        alpha=0.7,
                    )
                    # Also plot low, high (--) and MPV (|) lines vertically, if they fall within the histogram limits
                    # These vertical lines shouldn't have any legend entries and should be in the method_color
                    if (
                        histogram_limits[col][0]
                        <= low_point
                        <= histogram_limits[col][1]
                    ):
                        axs_processes[process][b, col].axvline(
                            low_point,
                            color=method_color[method],
                            linestyle="--",
                            alpha=0.7,
                        )
                    if (
                        histogram_limits[col][0]
                        <= high_point
                        <= histogram_limits[col][1]
                    ):
                        axs_processes[process][b, col].axvline(
                            high_point,
                            color=method_color[method],
                            linestyle="--",
                            alpha=0.7,
                        )
                    if (
                        histogram_limits[col][0]
                        <= mpv_point
                        <= histogram_limits[col][1]
                    ):
                        axs_processes[process][b, col].axvline(
                            mpv_point,
                            color=method_color[method],
                            linestyle=":",
                            alpha=0.7,
                        )
                    axs_processes[process][b, col].set_title(
                        f"[{bins_E[b]}, {bins_E[b+1]}] GeV"
                    )
                    # If it's the last column, put legend to the bottom of the plot
                    if (
                        col == len(histogram_limits) - 1
                        or col == len(histogram_limits) - 2
                    ):
                        axs_processes[process][b, col].legend(loc="lower right")
                    else:
                        axs_processes[process][b, col].legend()
                    axs_processes[process][b, col].grid(True)
        linestyle = "--"
        if method == "CaloJets_Durham":
            cp = f[process]["std68_all"]
        else:
            cp = f[process]["std68_all"]
        xs, ys = cp[2], cp[3]
        x_pts, y_pts = cp[4], cp[5]
        col, row = PROCESS_TO_ROW_COL[process]
        xs, ys, fit_params, fit_cov = get_func_fit(
            x_pts, y_pts, confusion_term="CaloJet" in method
        )
        ax[row, col].plot(
            xs,
            ys,
            label=f"{method_dict[method]} {print_params(fit_params)}",
            color=method_color[method],
            linestyle=linestyle,
        )
        ax[row, col].plot(x_pts, y_pts, "x", color=method_color[method], markersize=4)
        ax[row, col].set_ylabel("$\sigma_E / E_{true}$")
        ax[row, col].set_title(label)  # + rf"({print_noise_model(fit_params)})")
        if "CaloJet" in method:
            ax_fits = ax_fit_trials_calojets
        elif "IdealMatching" not in method:
            ax_fits = ax_fit_trials
        else:
            # skip
            continue

        ax_fits[row, col].plot(
            x_pts, y_pts, "x", markersize=4, label=method_dict[method]
        )
        # xs, ys, fit_params, fit_cov = get_func_fit(x_pts, y_pts, confusion_term=True)
        # ax_fits[row, col].plot(
        #    xs,
        #    ys,
        #    label=f"{print_params(fit_params)}",
        #    linestyle="--",
        # )
        xs, ys, fit_params, fit_cov = get_func_fit(
            x_pts,
            y_pts,
            confusion_term=True,
            error_summation=True,
            bounds_constant=[0.01, 0.03],
        )
        ax_fits[row, col].plot(
            xs,
            ys,
            label=f"{print_params(fit_params)}",
            linestyle="--",
        )
        # xs, ys, fit_params, fit_cov = get_func_fit(x_pts, y_pts, confusion_term=False)
        # ax_fits[row, col].plot(
        #    xs,
        #    ys,
        #    label=f"{print_params(fit_params)}",
        #    linestyle="--",
        # )
        xs, ys, fit_params, fit_cov = get_func_fit(
            x_pts,
            y_pts,
            confusion_term=False,
            error_summation=True,
            bounds_constant=[0.01, 0.03],
        )
        ax_fits[row, col].plot(
            xs,
            ys,
            label=f"{print_params(fit_params)}",
            linestyle="--",
        )
        xs, ys, fit_params, fit_cov = get_func_fit(
            x_pts,
            y_pts,
            confusion_term=True,
            bounds_constant=[0.005, 0.04],
            bounds_noise=[0.0, 0.6],
            bounds_confusion=[0.01, 0.6],
            min_E=0,
        )

        # ax_fits[row, col].plot(
        #    xs,
        #    ys,
        #    label=f" {print_params(fit_params)} (with bounds from 20GeV)",
        #    linestyle="--",
        # )
        ax_fits[row, col].plot(
            xs,
            ys,
            label=f"{print_params(fit_params)} (with bounds)",
            linestyle="--",
        )
        ax_fits[row, col].set_title(label)
        Higgs_x, Higgs_y = (
            f_mH[process]["x_vals_reco"],
            f_mH[process]["y_vals_reco"],
        )
        if LINE_STYLES[process] == ":":
            ax_mH[0, 0].step(
                Higgs_x,
                Higgs_y,
                where="post",
                label=f"{label} ({method_dict[method]})",
                color=PROCESS_COLORS[process],
                linestyle=method_linestyle[method],
            )
            ax_mH[0, 1].step(
                Higgs_x,
                Higgs_y,
                where="post",
                label=f"{label} ({method_dict[method]})",
                color=PROCESS_COLORS[process],
                linestyle=method_linestyle[method],
            )
            ax_mH[0, 0].set_title(r"H → Light-flavour jets")
        elif LINE_STYLES[process] == "-":
            ax_mH[1, 0].step(
                Higgs_x,
                Higgs_y,
                where="post",
                label=f"{label} ({method_dict[method]})",
                color=PROCESS_COLORS[process],
                linestyle=method_linestyle[method],
            )
            ax_mH[1, 1].step(
                Higgs_x,
                Higgs_y,
                where="post",
                label=f"{label} ({method_dict[method]})",
                color=PROCESS_COLORS[process],
                linestyle=method_linestyle[method],
            )
            ax_mH[1, 0].set_title(r"H → Light-flavour and b-jets")
for process in figs_processes:
    figs_processes[process].tight_layout()
    fig_process_path = os.path.join(
        outputDir, f"Detailed_JER_histograms_{process}_all.pdf"
    )
    print("Saving figure to", fig_process_path)
    figs_processes[process].savefig(fig_process_path)

for i in range(len(ax_mH)):
    for j in range(len(ax_mH[i])):
        ax_mH[i, j].set_xlabel("$m_H$ (reco.) [GeV]")
        ax_mH[i, j].set_xlabel("$m_H$ (reco.) [GeV]")
        ax_mH[i, j].grid()
    ax_mH[i, 1].set_xlim(100, 140)
    ax_mH[i, 0].legend(title="q ∈ {u, d, s}", fontsize=7.5, title_fontsize=8)
for i in range(len(ax)):
    for j in range(len(ax[i])):
        ax[i, j].legend(fontsize=6.5)
        ax[i, j].grid()
        ax[i, j].set_ylim([0, 0.35])
for ax in [ax_fit_trials, ax_fit_trials_calojets]:
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            ax[i, j].legend(fontsize=6.5)
            ax[i, j].grid()


fig.tight_layout()
annotate_matrix_plot_with_arrows(fig)
fig_path = os.path.join(outputDir, f"JER_comparison_PF_and_CaloJets.pdf")

print("Saving figure to", fig_path)
fig.savefig(fig_path)
fig_mH.tight_layout()
fig_mH_path = os.path.join(outputDir, f"Higgs_mass_comparison_PF_and_CaloJets.pdf")
print("Saving figure to", fig_mH_path)

fig_mH.savefig(fig_mH_path)
fig_E_mH_gluons.tight_layout()
fig_E_mH_gluons_path = os.path.join(
    outputDir, f"Higgs_mass_and_E_resolution_gluon_processes.pdf"
)
print("Saving figure to", fig_E_mH_gluons_path)
fig_E_mH_gluons.savefig(fig_E_mH_gluons_path)

fig_fit_trials_calojets.tight_layout()
annotate_matrix_plot_with_arrows(fig_fit_trials_calojets)
path_cj = os.path.join(outputDir, f"JER_fitting_CaloJets.pdf")
fig_fit_trials_calojets.savefig(path_cj)

# Similar for fig_fit_trials
fig_fit_trials.tight_layout()
annotate_matrix_plot_with_arrows(fig_fit_trials)
path_pf = os.path.join(outputDir, f"JER_fitting_PFJets.pdf")
fig_fit_trials.savefig(path_pf)

if not args.AK_comparison:
    for prefix in ["neutral", "charged", "photons"]:
        ###########################################
        # Resolution plots, calo jets vs. neutral part of PF jets
        fig, ax = plt.subplots(5, 3, figsize=(15, 12))
        for method in [
            "PF_Durham",
            "PF_Durham_IdealMatching",
            "CaloJets_Durham",
        ]:
            method_color = {
                "PF_Durham": "blue",
                "CaloJets_Durham": "green",
                "PF_Durham_IdealMatching": "orange",
            }
            method_linestyle = {  # For the $m_H$ plots
                "PF_Durham": "-",
                "CaloJets_Durham": "--",
                "PF_Durham_IdealMatching": ":",
            }
            method_text = {
                "PF_Durham": "PF jets, {}".format(prefix),
                "CaloJets_Durham": "Calo jets",
                "PF_Durham_IdealMatching": "PF jets, ideal matching",
            }
            f = pickle.load(
                open(
                    os.path.join(
                        inputDir,
                        f"{method}/plots_resolution/energy_fit_params_per_process.pkl",
                    ),
                    "rb",
                )
            )
            f_mH = pickle.load(
                open(
                    os.path.join(
                        inputDir, f"{method}/plots_mass/Higgs_mass_histograms_data.pkl"
                    ),
                    "rb",
                )
            )
            for i, process in enumerate(PROCESS_COLORS.keys()):
                if process not in f:
                    continue
                color = PROCESS_COLORS[process]
                label = HUMAN_READABLE_PROCESS_NAMES[process]
                linestyle = "--"
                if method == "CaloJets_Durham":
                    cp = f[process]["std68_all"]
                else:
                    cp = f[process]["std68_" + prefix]
                # xs, ys = cp[2], cp[3]
                # estimate parameters here from x_pts and y_pts
                x_pts, y_pts = cp[4], cp[5]
                xs, ys, fit_params, fit_cov = get_func_fit(x_pts, y_pts)
                col, row = PROCESS_TO_ROW_COL[process]
                ax[row, col].plot(
                    xs,
                    ys,
                    label=f"{method_text[method]} {print_params(fit_params)}",
                    color=method_color[method],
                    linestyle=linestyle,
                )
                ax[row, col].plot(
                    x_pts, y_pts, "x", color=method_color[method], markersize=4
                )
                ax[row, col].set_ylabel("$\sigma_E / E_{true}$")
                ax[row, col].set_title(label)
                ax[row, col].set_xlabel("$E_{true}$ [GeV]")
                ax[row, col].set_ylim([0, 0.35])
                Higgs_x, Higgs_y = (
                    f_mH[process]["x_vals_reco"],
                    f_mH[process]["y_vals_reco"],
                )
        for i in range(len(ax)):
            for j in range(len(ax[i])):
                ax[i, j].legend(fontsize=6.5)
                ax[i, j].grid()
        fig.tight_layout()
        annotate_matrix_plot_with_arrows(fig)
        fig_path = os.path.join(
            outputDir, f"JER_comparison_{prefix}_PF_vs_CaloJets.pdf"
        )
        print("Saving figure to", fig_path)
        fig.savefig(fig_path)
