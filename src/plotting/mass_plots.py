# For each root file in the direct inputDir, open the root histogram and read the 'h_fancy' histogram. the legend entry should be the root file name. plot it on the same mpl canvas and please normalize it to 1!
import os
import ROOT
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from src.process_config import (
    PROCESS_COLORS,
    HUMAN_READABLE_PROCESS_NAMES,
    LINE_STYLES,
    NUMBER_OF_JETS,
)
import pickle
import argparse


matplotlib.rcParams.update(
    {
        #'font.sans-serif': "Arial",
        "font.family": "sans-serif",  # Ensure Matplotlib uses the sans-serif family
        # "mathtext.fontset": "stix", # serif math, similar to LaTeX Times
        # "mathtext.default": "it",   # math variables italic by default
        "font.size": 12,
    }
)


parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", type=str, required=True)

args = parser.parse_args()

inputDir = args.inputDir
outputDir = os.path.join(inputDir, "plots_mass")
os.makedirs(outputDir, exist_ok=True)

root_files = [f for f in os.listdir(inputDir) if f.endswith(".root")]


fig, ax = plt.subplots(3, 1, figsize=(6, 9))  #
for fname in sorted(root_files):
    file_path = os.path.join(inputDir, fname)
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue
    hist = f.Get("h_jet_E_reco_over_E_true")
    if not hist:
        print(f"No 'h_jet_E_reco_over_E_true' histogram in {fname}")
        f.Close()
        continue
    # Convert histogram to numpy arrays
    n_bins = hist.GetNbinsX()
    x_vals = np.array([hist.GetBinCenter(i) for i in range(1, n_bins + 1)])
    y_vals = np.array([hist.GetBinContent(i) for i in range(1, n_bins + 1)])
    # Normalize
    integral = np.sum(y_vals)
    print("Integral of histogram in {}: {}".format(fname, integral))
    if integral > 0:
        y_vals = y_vals / integral
    else:
        print(f"Warning: {fname} histogram integral = 0")
    # Plot
    label = os.path.splitext(fname)[0]
    # plt.plot(x_vals, y_vals, "x", label=HUMAN_READABLE_PROCESS_NAMES[label], color=PROCESS_COLORS[label])
    plt.plot(
        x_vals,
        y_vals,
        linestyle=LINE_STYLES[label],
        color=PROCESS_COLORS[label],
        label=HUMAN_READABLE_PROCESS_NAMES[label],
    )
    f.Close()
plt.xlabel("$E_{reco} / E_{true}")
plt.xlim([0.8, 1.2])
plt.ylabel("Normalized Entries")
# plt.title()
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(inputDir, "norm_E_over_true_overlaid.pdf"))

# Make the same plot but with a thin line, only for nu nu q q and with xlim from -0.6 to +0.6
fig, ax = plt.subplots(3, 1, figsize=(6, 10.5))
for fname in root_files:
    file_path = os.path.join(inputDir, fname)
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue
    hist = f.Get("h_jet_E_reco_over_E_true")
    if not hist:
        print(f"No 'h_jet_E_reco_over_E_true' histogram in {fname}")
        f.Close()
        continue
    # Convert histogram to numpy arrays
    n_bins = hist.GetNbinsX()
    x_vals = np.array([hist.GetBinCenter(i) for i in range(1, n_bins + 1)])
    y_vals = np.array([hist.GetBinContent(i) for i in range(1, n_bins + 1)])
    # Normalize
    integral = np.sum(y_vals)
    print("Integral of histogram in {}: {}".format(fname, integral))
    if integral > 0:
        y_vals = y_vals / integral
    else:
        print(f"Warning: {fname} histogram integral = 0")
    # Plot
    label = os.path.splitext(fname)[0]
    ax[0].plot(
        x_vals,
        y_vals,
        label=HUMAN_READABLE_PROCESS_NAMES[label],
        linewidth=1,
        linestyle=LINE_STYLES[label],
        color=PROCESS_COLORS[label],
    )
    ax[1].plot(
        x_vals,
        y_vals,
        label=HUMAN_READABLE_PROCESS_NAMES[label],
        linewidth=1,
        linestyle=LINE_STYLES[label],
        color=PROCESS_COLORS[label],
    )
    ax[2].plot(
        x_vals,
        y_vals,
        label=HUMAN_READABLE_PROCESS_NAMES[label],
        linewidth=1,
        linestyle=LINE_STYLES[label],
        color=PROCESS_COLORS[label],
    )
    # Now make a Gaussian fit from x_vals 0.8 to 1.2. Plot it in that range as well and put sigma and mean in the legend
    # You can initialize parameters with stdev and mean
    # Do it here:
    # Fit a Gaussian to the 0.8–1.2 range and plot it
    mask = (x_vals >= 0.8) & (x_vals <= 1.2)
    x_fit = x_vals[mask]
    y_fit = y_vals[mask]
    if x_fit.size >= 3 and np.sum(y_fit) > 0:
        # Moment-based initial guesses (use y as weights)
        w = y_fit
        mu0 = np.average(x_fit, weights=w)
        var0 = np.average((x_fit - mu0) ** 2, weights=w)
        sigma0 = np.sqrt(max(var0, 1e-12))
        A0 = y_fit.max()

        def gauss(x, A, mu, sigma):
            # clip sigma to avoid division by zero during fitting/plotting
            s = np.clip(sigma, 1e-12, None)
            return A * np.exp(-0.5 * ((x - mu) / s) ** 2)

        # Try a nonlinear least-squares fit; fall back to moment estimates if it fails
        try:
            from scipy.optimize import curve_fit

            popt, _ = curve_fit(
                gauss,
                x_fit,
                y_fit,
                p0=[A0, mu0, sigma0],
                bounds=([0.0, 0.8, 1e-6], [np.inf, 1.2, np.inf]),
                maxfev=10000,
            )
            A, mu, sigma = popt
        except Exception as _e:
            # Fallback: use moment estimates without optimization
            A, mu, sigma = A0, mu0, sigma0
        # Plot the fitted Gaussian over the fit window
        x_dense = np.linspace(0.9, 1.1, 200)
        ax[0].plot(
            x_dense,
            gauss(x_dense, A, mu, sigma),
            linestyle="--",
            linewidth=1,
        )
        ax[1].plot(
            x_dense,
            gauss(x_dense, A, mu, sigma),
            linestyle="--",
            linewidth=1,
        )
        ax[2].plot(
            x_dense,
            gauss(x_dense, A, mu, sigma),
            linestyle="--",
            linewidth=1,
        )
    else:
        print(f"Not enough points in [0.8, 1.2] for {fname} to fit.")
    for i in range(3):
        ax[i].set_xlabel("$E_{reco}$ / $E_{true}$")
        ax[i].set_ylabel("Jets (normalized)")
    # ax[0].set_title()
    ax[0].legend()
    ax[0].grid(True)
    ax[1].grid(True)
    ax[2].grid(True)
    ax[1].set_xlim([0.85, 1.15])
    ax[2].set_xlim([0.95, 1.05])
    f.Close()

fig.tight_layout()
fig.savefig(os.path.join(outputDir, "norm_E_over_true_overlaid_v2.pdf"))

# There are two histograms: h_genjet_all_energies and h_genjet_matched_energies. Make a plot with the ratio (so basically efficiency) of matched over all vs energy
plt.figure(figsize=(6, 6))
for fname in sorted(root_files):
    file_path = os.path.join(inputDir, fname)
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue
    hist_all = f.Get("h_genjet_all_energies")
    hist_matched = f.Get("h_genjet_matched_energies")
    if not hist_all or not hist_matched:
        print(f"No required histograms in {fname}")
        f.Close()
        continue
    # Convert histograms to numpy arrays
    n_bins = hist_all.GetNbinsX()
    x_vals = np.array([hist_all.GetBinCenter(i) for i in range(1, n_bins + 1)])
    # Remove the xvals larger than 175
    filt = (
        np.array([hist_all.GetBinContent(i) for i in range(1, n_bins + 1)]) > 1000
    )  # Cut out the low statistics bins
    filt = filt & (x_vals <= 100)
    x_vals = x_vals[filt]
    y_all = np.array([hist_all.GetBinContent(i) for i in range(1, n_bins + 1)])[filt]
    y_matched = np.array([hist_matched.GetBinContent(i) for i in range(1, n_bins + 1)])[
        filt
    ]
    # Calculate ratio
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.true_divide(y_matched, y_all)
        ratio[~np.isfinite(ratio)] = 0  # set inf and NaN to 0
    # PlotR
    label = os.path.splitext(fname)[0]
    plt.plot(
        x_vals,
        ratio,
        LINE_STYLES[label],
        label=HUMAN_READABLE_PROCESS_NAMES[label],
        color=PROCESS_COLORS[label],
    )
    plt.plot(x_vals, ratio, "x", color=PROCESS_COLORS[label])
    f.Close()

plt.xlabel("$E_{true}$ [GeV]")
plt.ylabel("Matching Efficiency (Matched/ All)")
# plt.ylim([0.80, 1.001])
plt.title("Jet Matching Efficiency vs. Energy")
plt.legend(title="q ∈ {u, d, s}", fontsize=10, title_fontsize=8)
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(outputDir, "matching_efficiency_vs_energy.pdf"))
plt.clf()


## Produce the Higgs mass plots similar to the ones in the fccanalysis plots file but easier to manipulate
fig, ax = plt.subplots(len(root_files), 2, figsize=(8, 2.7 * len(root_files)))
fig_mH_all, ax_mH_all = plt.subplots(
    1, 2, figsize=(8, 6)
)  # plot mH reco of all root files on the same plot, left side fully and right side zoomed into the peak

fig_mH_njets, ax_mH_njets = plt.subplots(
    3, 2, figsize=(7, 9)
)  # Plot separately: processes with 2 jets, 4 jets, 6 jets in the final state (each row would have a full plot of like 50-175 GeV and the second one zoomed in 105-145 GeV)
fig_mH_perprocess, ax_mH_perprocess = plt.subplots(
    2, 2, figsize=(7, 9)
)  # First plot: 2, 4, 6 jets for the full lines (heavy flavour); second plot: 2,4,6 jets for dotted lines (light flavour)

figlog, axlog = plt.subplots(len(root_files), 2, figsize=(8, 2.7 * len(root_files)))

if len(root_files) == 1:
    ax = np.array([ax])
    axlog = np.array([axlog])

process_to_mH_hist_plots = {}

# Higgs mass histogram
for i, fname in enumerate(sorted(root_files)):
    file_path = os.path.join(inputDir, fname)
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue
    hist_gen = f.Get("h_mH_gen")
    hist_reco = f.Get("h_mH_reco")
    hist_gt = f.Get(
        "h_mH_stable_gt_particles"
    )  # for the resolution with ideal clustering
    hist_gt_recomatched = f.Get("h_mH_reco_particles_matched")
    if not hist_gen or not hist_reco:
        print(f"No 'h_mH_gen'/'h_mH_reco' histogram in {fname}")
        f.Close()
        continue
    n_bins = hist_gen.GetNbinsX()
    x_vals_gen = np.array([hist_gen.GetBinCenter(i) for i in range(1, n_bins + 1)])
    y_vals_gen = np.array([hist_gen.GetBinContent(i) for i in range(1, n_bins + 1)])

    nbins_gt = hist_gt.GetNbinsX()
    x_vals_gt = np.array([hist_gt.GetBinCenter(i) for i in range(1, nbins_gt + 1)])
    y_vals_gt = np.array([hist_gt.GetBinContent(i) for i in range(1, nbins_gt + 1)])

    nbins_gt_reco_matched = hist_gt_recomatched.GetNbinsX()
    x_vals_gt_recomatched = np.array(
        [
            hist_gt_recomatched.GetBinCenter(i)
            for i in range(1, nbins_gt_reco_matched + 1)
        ]
    )
    y_vals_gt_recomatched = np.array(
        [
            hist_gt_recomatched.GetBinContent(i)
            for i in range(1, nbins_gt_reco_matched + 1)
        ]
    )
    step_size_gt_recomatched = x_vals_gt_recomatched[1] - x_vals_gt_recomatched[0]
    integral_gt_recomatched = np.sum(y_vals_gt_recomatched)
    print(
        "Integral of reco-GT matched histogram in {}: {}".format(
            fname, integral_gt_recomatched
        )
    )
    if integral_gt_recomatched > 0:
        y_vals_gt_recomatched = (
            y_vals_gt_recomatched / integral_gt_recomatched / step_size_gt_recomatched
        )
    else:
        print(f"Warning: {fname} reco-GT matched histogram integral = 0")
    # Normalize
    integral = np.sum(y_vals_gen)
    step_size_gen = x_vals_gen[1] - x_vals_gen[0]
    print("Integral of histogram in {}: {}".format(fname, integral))
    if integral > 0:
        y_vals_gen = y_vals_gen / integral / step_size_gen
    else:
        print(f"Warning: {fname} histogram integral = 0")
    print("Sum of y vals now", np.sum(y_vals_gen))
    integral_gt = np.sum(y_vals_gt)
    step_size_gt = x_vals_gt[1] - x_vals_gt[0]
    print("Integral of GT histogram in {}: {}".format(fname, integral_gt))
    if integral_gt > 0:
        y_vals_gt = y_vals_gt / integral_gt / step_size_gt
    else:
        print(f"Warning: {fname} GT histogram integral = 0")
    # Plot
    label = os.path.splitext(fname)[0]
    # Plot the histograms with bins for reco (full lines) and gen(dashed lines)
    n_bins = hist_reco.GetNbinsX()
    x_vals_reco = np.array([hist_reco.GetBinCenter(i) for i in range(1, n_bins + 1)])
    y_vals_reco = np.array([hist_reco.GetBinContent(i) for i in range(1, n_bins + 1)])
    integral_reco = np.sum(y_vals_reco)
    step_size_reco = x_vals_reco[1] - x_vals_reco[0]
    print("Integral of histogram in {}: {}".format(fname, integral_reco))
    if integral_reco > 0:
        y_vals_reco = y_vals_reco / integral_reco
        y_vals_reco = y_vals_reco / step_size_reco  # Normalize to bin width
    else:
        print(f"Warning: {fname} histogram integral = 0")
    # ax_mH_all[0].step(x_vals_reco, y_vals_reco, where='mid', label=label, linestyle='solid')
    # ax_mH_all[1].step(x_vals_reco, y_vals_reco, where='mid', label=label, linestyle='solid')

    ax_mH_all[0].step(
        x_vals_reco,
        y_vals_reco,
        where="mid",
        color=PROCESS_COLORS[label],
        linestyle=LINE_STYLES[label],
        label=HUMAN_READABLE_PROCESS_NAMES[label],
    )
    ax_mH_all[1].step(
        x_vals_reco,
        y_vals_reco,
        where="mid",
        color=PROCESS_COLORS[label],
        linestyle=LINE_STYLES[label],
        label=HUMAN_READABLE_PROCESS_NAMES[label],
    )

    process_to_mH_hist_plots[label] = {
        "x_vals_reco": x_vals_reco,
        "y_vals_reco": y_vals_reco,
        "x_vals_gen": x_vals_gen,
        "y_vals_gen": y_vals_gen,
        "x_vals_gt": x_vals_gt,
        "y_vals_gt": y_vals_gt,
        "x_vals_gt_recomatched": x_vals_gt_recomatched,
        "y_vals_gt_recomatched": y_vals_gt_recomatched,
    }

    # Plot the step onto ax_mH_NJets based on NUMBER_OF_JETS
    njets = NUMBER_OF_JETS.get(label, None)
    if njets in [2, 4, 6]:
        row = {2: 0, 4: 1, 6: 2}[njets]
        ax_mH_njets[row, 0].step(
            x_vals_reco,
            y_vals_reco,
            where="mid",
            label=HUMAN_READABLE_PROCESS_NAMES[label],
        )
        ax_mH_njets[row, 1].step(
            x_vals_reco,
            y_vals_reco,
            where="mid",
            label=HUMAN_READABLE_PROCESS_NAMES[label],
        )
        ax_mH_njets[row, 1].set_xlim([80, 150])
        ax_mH_njets[row, 0].set_xlabel("$m_H$ [GeV]")
        ax_mH_njets[row, 1].set_xlabel("$m_H$ [GeV]")
        ax_mH_njets[row, 0].set_ylabel("Normalized Events")
        ax_mH_njets[row, 0].set_title(f"{njets} jets")
    # Plot the step onto ax_mH_perprocess based on line style (full line = heavy flavour, dotted = light flavour)
    if LINE_STYLES.get(label, "") == "-":  # full line
        row = 0
    elif LINE_STYLES.get(label, "") == ":":  # dotted line
        row = 1
    else:
        row = None
    if row is not None:
        ax_mH_perprocess[row, 0].step(
            x_vals_reco,
            y_vals_reco,
            where="mid",
            color=PROCESS_COLORS[label],
            linestyle="-",
            label=HUMAN_READABLE_PROCESS_NAMES[label],
        )
        ax_mH_perprocess[row, 1].step(
            x_vals_reco,
            y_vals_reco,
            where="mid",
            color=PROCESS_COLORS[label],
            linestyle="-",
            label=HUMAN_READABLE_PROCESS_NAMES[label],
        )
        ax_mH_perprocess[row, 1].set_xlim([80, 150])
        ax_mH_perprocess[row, 0].set_xlabel("$m_H$ [GeV]")
        ax_mH_perprocess[row, 1].set_xlabel("$m_H$ [GeV]")
        ax_mH_perprocess[row, 0].set_ylabel("Normalized Events")
        if row == 0:
            ax_mH_perprocess[row, 0].set_title("Containing b-jets")
        else:
            ax_mH_perprocess[row, 0].set_title("Containing only light-flavour jets")
    # ax.plot(x_vals_reco, y_vals_reco, label=label + " (reco)", linestyle='solid')
    # ax.plot(x_vals_gen, y_vals_gen, label=label + " (gen)", linestyle='dashed')
    for k in range(2):
        ax[i, k].step(x_vals_reco, y_vals_reco, where="mid", color="blue", label="reco")
        ax[i, k].step(x_vals_gen, y_vals_gen, where="mid", color="orange", label="gen")
        # ax[i].step(x_vals_gt, y_vals_gt, where='mid', color="green", label="GT", linestyle='dotted')
        axlog[i, k].step(
            x_vals_reco, y_vals_reco, where="mid", color="blue", label="reco"
        )
        axlog[i, k].step(
            x_vals_gen, y_vals_gen, where="mid", color="orange", label="gen"
        )
        axlog[i, k].step(x_vals_gt, y_vals_gt, where="mid", color="green", label="GT")
        axlog[i, k].step(
            x_vals_gt_recomatched,
            y_vals_gt_recomatched,
            where="mid",
            color="red",
            label="reco-GT matched",
        )
        axlog[i, k].legend()
        axlog[i, k].set_ylim([1e-3, 1])
        ax[i, k].legend()
        axlog[i, k].set_title(HUMAN_READABLE_PROCESS_NAMES[label])
        axlog[i, k].set_yscale("log")
        axlog[i, k].set_xlabel("$m_H$ [GeV]")
        axlog[i, k].set_ylabel("Events (norm.)")
        ax[i, k].set_title(HUMAN_READABLE_PROCESS_NAMES[label])
    ax[i, 1].set_xlim([80, 140])
    axlog[i, 1].set_xlim([80, 140])

p = os.path.join(outputDir, "Higgs_mass_reco_vs_gen.pdf")
plog = os.path.join(outputDir, "log_Higgs_mass_reco_vs_gen.pdf")
phiggs = os.path.join(outputDir, "Higgs_mass_reco_overlaid_mH_reco_normalized.pdf")
path_higgs_separate_by_process_type = os.path.join(
    outputDir, "Higgs_mass_plots_sorted_per_process_type.pdf"
)
path_higgs_separate_by_njets = os.path.join(
    outputDir, "Higgs_mass_plots_sorted_per_N_jets.pdf"
)
path_Higgs_pkl = os.path.join(outputDir, "Higgs_mass_histograms_data.pkl")

pickle.dump(process_to_mH_hist_plots, open(path_Higgs_pkl, "wb"))

for i in range(len(ax_mH_perprocess)):
    ax_mH_perprocess[i, 0].legend(title="q ∈ {u, d, s}", fontsize=11, title_fontsize=9)
    ax_mH_perprocess[i, 0].grid()
    ax_mH_perprocess[i, 1].grid()

for i in range(len(ax_mH_njets)):
    ax_mH_njets[i, 0].legend(title="q ∈ {u, d, s}", fontsize=11, title_fontsize=9)
    ax_mH_njets[i, 0].grid()
    ax_mH_njets[i, 1].grid()

ax_mH_all[0].set_xlabel("$m_H$ [GeV]")
ax_mH_all[1].set_xlabel("$m_H$ [GeV]")
ax_mH_all[0].set_ylabel("Normalized Events")
ax_mH_all[1].set_ylabel("Normalized Events")
ax_mH_all[1].set_xlim([115, 135])
ax_mH_all[0].grid()
ax_mH_all[1].grid()
ax_mH_all[0].legend(title="q ∈ {u, d, s}", fontsize=12, title_fontsize=10)
# ax_mH_all[1].legend()

fig_mH_all.tight_layout()
fig.tight_layout()
fig.savefig(p)
figlog.tight_layout()
figlog.savefig(plog)
fig_mH_all.savefig(phiggs)
print("Saving to", p, plog, phiggs)

fig_mH_perprocess.tight_layout()
fig_mH_njets.tight_layout()
fig_mH_perprocess.savefig(path_higgs_separate_by_process_type)
fig_mH_njets.savefig(path_higgs_separate_by_njets)
print("Saving to", path_higgs_separate_by_process_type, path_higgs_separate_by_njets)
plt.clf()


# Make a plot of h_mH_all_stable_part
fig2, ax2 = plt.subplots(len(root_files), 1, figsize=(8, 2.7 * len(root_files)))
if len(root_files) == 1:
    ax2 = np.array([ax2])
for i, fname in enumerate(root_files):
    file_path = os.path.join(inputDir, fname)
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue
    hist_gen = f.Get("h_mH_all_stable_part")
    if not hist_gen:
        print(f"No 'h_mH_all_stable_part' histogram in {fname}")
        f.Close()
        continue
    n_bins = hist_gen.GetNbinsX()
    x_vals_gen = np.array([hist_gen.GetBinCenter(i) for i in range(1, n_bins + 1)])
    y_vals_gen = np.array([hist_gen.GetBinContent(i) for i in range(1, n_bins + 1)])
    # Normalize
    integral = np.sum(y_vals_gen)
    step_size_gen = x_vals_gen[1] - x_vals_gen[0]
    print("Integral of histogram in {}: {}".format(fname, integral))
    if integral > 0:
        y_vals_gen = y_vals_gen / integral / step_size_gen
    else:
        print(f"Warning: {fname} histogram integral = 0")
    # Plot
    label = os.path.splitext(fname)[0]
    ax2[i].step(x_vals_gen, y_vals_gen, where="mid", label=label)
    ax2[i].set_title(label)
    ax2[i].set_xlabel("Inv Mass all gen particles [GeV]")
    ax2[i].set_ylabel("Normalized Entries / GeV")
    ax2[i].legend()
    ax2[i].grid()
    f.Close()
fig2.tight_layout()
fig2.savefig(os.path.join(outputDir, "inv_mass_all_gen_particles_normalized.pdf"))
plt.clf()
# plt.show()

# Make a histogram of h_E_reco_over_true_Charged. First, zoomed in from 0.9 to 1.1, then, log y from 0 to 2 (full range of the histogram)
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
for fname in root_files:
    file_path = os.path.join(inputDir, fname)
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue
    hist = f.Get("h_E_reco_over_true_Charged")
    if not hist:
        print(f"No 'h_E_reco_over_true_Charged' histogram in {fname}")
        f.Close()
        continue
    n_bins = hist.GetNbinsX()
    x_vals = np.array([hist.GetBinCenter(i) for i in range(1, n_bins + 1)])
    y_vals = np.array([hist.GetBinContent(i) for i in range(1, n_bins + 1)])
    # Normalize
    integral = np.sum(y_vals)
    step_size = x_vals[1] - x_vals[0]
    print("Integral of histogram in {}: {}".format(fname, integral))
    if integral > 0:
        y_vals = y_vals / integral / step_size
    else:
        print(f"Warning: {fname} histogram integral = 0")
    label = os.path.splitext(fname)[0]
    ax[0].step(x_vals, y_vals, where="mid", label=label)
    ax[1].step(x_vals, y_vals, where="mid", label=label)
    f.Close()

ax[0].set_xlabel("E_reco / E_true (Charged)")
ax[0].set_xlim([0.9, 1.1])
ax[0].legend()
ax[1].set_xlabel("E_reco / E_true (Charged)")
ax[1].set_yscale("log")
ax[1].legend()
ax[0].set_ylabel("Normalized Entries / bin width")
fig.tight_layout()
fig.savefig(os.path.join(outputDir, "E_reco_over_true_Charged.pdf"))
plt.clf()

# Histogram of ratio_jet_energies_fancy_Neutral_part and ratio_jet_energies_fancy_Charged_part and ratio_jet_energies_fancy
# so basically neutral, charged and all particles - make one histogram for each root file (with three lines, Neutral, Charged, and All) and stack them vertically
binsE = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]

fig, ax = plt.subplots(
    len(root_files), len(binsE), figsize=(4 * len(binsE), 4 * len(root_files))
)

if len(root_files) == 1:
    ax = np.array([ax])

for i, fname in enumerate(root_files):
    print(
        "===> Plotting binned jet-level Ereco/Etrue neutral/charged/all for file", fname
    )
    for j in range(len(binsE) - 1):
        file_path = os.path.join(inputDir, fname)
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie():
            print(f"Could not open {fname}")
            continue
        hist_neutral = f.Get(
            "binned_E_Neutral_reco_over_true_{}_{}".format(binsE[j], binsE[j + 1])
        )
        hist_charged = f.Get(
            "binned_E_Charged_reco_over_true_{}_{}".format(binsE[j], binsE[j + 1])
        )
        hist_all = f.Get("binned_E_reco_over_true_{}_{}".format(binsE[j], binsE[j + 1]))
        n_bins = hist_all.GetNbinsX()
        x_vals = np.array([hist_all.GetBinCenter(i) for i in range(1, n_bins + 1)])
        if hist_neutral:
            y_vals_neutral = np.array(
                [hist_neutral.GetBinContent(i) for i in range(1, n_bins + 1)]
            )
            y_vals_charged = np.array(
                [hist_charged.GetBinContent(i) for i in range(1, n_bins + 1)]
            )
            integral_neutral = np.sum(y_vals_neutral)
            integral_charged = np.sum(y_vals_charged)
        y_vals_all = np.array([hist_all.GetBinContent(i) for i in range(1, n_bins + 1)])
        # Normalize
        integral_all = np.sum(y_vals_all)
        step_size = x_vals[1] - x_vals[0]
        if hist_neutral:
            print(
                "Integrals in {}: Neutral {}, Charged {}, All {}".format(
                    fname, integral_neutral, integral_charged, integral_all
                )
            )
            if integral_neutral > 0:
                y_vals_neutral = y_vals_neutral / integral_neutral / step_size
            else:
                print(f"Warning: {fname} neutral histogram integral = 0")
            if integral_charged > 0:
                y_vals_charged = y_vals_charged / integral_charged / step_size
            else:
                print(f"Warning: {fname} charged histogram integral = 0")
        if integral_all > 0:
            y_vals_all = y_vals_all / integral_all / step_size
        else:
            print(f"Warning: {fname} all histogram integral = 0")
        label = os.path.splitext(fname)[0]
        if hist_neutral:
            ax[i, j].step(
                x_vals, y_vals_neutral, where="mid", label="Neutral", color="red"
            )
            ax[i, j].step(
                x_vals, y_vals_charged, where="mid", label="Charged", color="blue"
            )
        ax[i, j].step(x_vals, y_vals_all, where="mid", label="All", color="green")
        ax[i, j].set_xlim([0.9, 1.1])
        ax[i, j].set_title(label)
        ax[i, j].set_xlabel("E_reco / E_true")
        ax[i, j].set_ylabel("Normalized Entries / bin width")
        ax[i, j].legend()
        ax[i, j].grid()
        label = os.path.splitext(fname)[0]
        ax[i, j].set_title(
            "{} [{},{}]GeV".format(
                HUMAN_READABLE_PROCESS_NAMES[label], binsE[j], binsE[j + 1]
            )
        )
        f.Close()

fig.tight_layout()
fig.savefig(os.path.join(outputDir, "ratio_jet_energies_fancy_Neutral_Charged_All.pdf"))

# Now, convert all ax[i, j] into log y scale and then save this as the same file + "_log" in the end...
for i in range(len(root_files)):
    for j in range(len(binsE)):
        ax[i, j].set_xlim([0.9, 1.1])
        ax[i, j].set_yscale("log")

fig.tight_layout()
fig.savefig(
    os.path.join(outputDir, "ratio_jet_energies_fancy_Neutral_Charged_All_LOGSCALE.pdf")
)
plt.clf()
