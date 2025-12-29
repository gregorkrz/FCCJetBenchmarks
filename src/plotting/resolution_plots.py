import numpy as np
import ROOT
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from copy import copy
import argparse
import os
from src.process_config import PROCESS_COLORS, HUMAN_READABLE_PROCESS_NAMES, LINE_STYLES
import pickle
import matplotlib

matplotlib.rcParams.update(
    {
        #'font.sans-serif': "Arial",
        "font.family": "sans-serif",  # Ensure Matplotlib uses the sans-serif family
        # "mathtext.fontset": "stix",           # serif math, similar to LaTeX Times
        # "mathtext.default": "it",             # math variables italic by default
        "font.size": 12,
    }
)

parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", type=str, required=True)
parser.add_argument("--angles-only", action="store_true", help="Only plot angular resolutions")

args = parser.parse_args()

inputDir = args.inputDir

outputDir = os.path.join(args.inputDir, "plots_resolution")
# Make dir if it doesn't exist

os.makedirs(outputDir, exist_ok=True)
print("Saving to directory:", outputDir)


def print_params(popt):
    if len(popt) == 2:
        return f"A={round(popt[0], 2)} C={round(popt[1], 2)}"
    return f"A={round(popt[0], 2)} B={round(popt[2], 2)} C={round(popt[1], 2)}"


def point_format(number):
    return str(number).replace(".", "p")


def neg_format(number):
    # put n5 for -5
    if number < 0:
        return point_format("n{}".format(abs(number)))
    else:
        return point_format(number)


processList = {}

for file in os.listdir(inputDir):
    if file.endswith(".root"):
        proc_name = file.replace(".root", "")
        processList[proc_name] = {"fraction": 1}

processList = {
    "p8_ee_ZH_vvqq_ecm240": {"fraction": 1},
}

########################################################################################################

binsE = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
bins_eta = [-5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 5]
bins_costheta = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]


# Define the Double-Sided Crystal Ball function
def double_crystal_ball(x, mu, sigma, alphaL, nL, alphaR, nR, norm):
    """Double-sided Crystal Ball PDF."""
    t = (x - mu) / sigma
    result = np.zeros_like(t)

    # Left side
    maskL = t < -alphaL
    result[maskL] = norm * (
        (nL / abs(alphaL)) ** nL
        * np.exp(-0.5 * alphaL**2)
        / (nL / abs(alphaL) - abs(alphaL) - t[maskL]) ** nL
    )

    # Gaussian core
    maskC = (~maskL) & (t < alphaR)
    result[maskC] = norm * np.exp(-0.5 * t[maskC] ** 2)

    # Right side
    maskR = t >= alphaR
    result[maskR] = norm * (
        (nR / abs(alphaR)) ** nR
        * np.exp(-0.5 * alphaR**2)
        / (nR / abs(alphaR) - abs(alphaR) + t[maskR]) ** nR
    )
    return result


def get_result_for_process(
    procname,
    bins=binsE,
    suffix="",
    sigma_method="std68",
    root_histogram_prefix="binned_E_reco_over_true_",
    wmin=0.7,
    wmax=1.2,
    divide_by_MPV=True
):
    # Sigma methods: std68, RMS, interquantile_range
    field_names = []
    f = ROOT.TFile.Open(os.path.join(inputDir, "{}.root".format(procname)))
    fig_hist, ax_hist = plt.subplots(
        3, 1, figsize=(8, 8.5)
    )  # stack two plots vertically
    eps = 0.0005
    if not divide_by_MPV:
        eps = 0.001 # try a bit bigger epsilon for angle resolution
    # else:
    #    raise Exception # TEMPORARILY

    def get_std68(
        theHist, bin_edges, percentage=0.683, epsilon=eps, wmin=0.7, wmax=1.2
    ):
        # theHist, bin_edges = np.histogram(data_for_hist, bins=bins, density=True)
        if wmin > 0:
            theHist[0] = 0.0  # for the energy histograms
        s = np.sum(theHist * np.diff(bin_edges))
        if s != 0:
            theHist /= s  # normalize the histogram to 1
        weight = 0.0
        #print("Bin edges:", bin_edges[:10], bin_edges[-10:], "len", len(bin_edges))
        #print("Histogram:", theHist[:10], theHist[-10:],"max",  max(theHist), "len", len(theHist))
        points = []
        sums = []
        am = np.argmax(theHist)
        MPV = 0.5 * (bin_edges[am] + bin_edges[am + 1])
        # Fill list of bin centers and the integral up to those point
        for i in range(len(bin_edges) - 1):
            weight += theHist[i] * (bin_edges[i + 1] - bin_edges[i])
            sums.append(weight)
            if bin_edges[i + 1] < wmin or bin_edges[i] > wmax:
                continue
            points.append([(bin_edges[i + 1] + bin_edges[i]) / 2, weight])
        low = wmin
        high = wmax
        width = 100
        for i in range(len(points)):
            for j in range(i, len(points)):
                wy = points[j][1] - points[i][1]
                if abs(wy - percentage) < epsilon:
                    wx = points[j][0] - points[i][0]
                    if wx < width:
                        low = points[i][0]
                        high = points[j][0]
                        width = wx
        if low == wmin and high == wmax:
            # Didn't fit well, try mean and stdev
            # Compute the stdev from the histogram
            print("Fitting didn't work")
            mean = np.sum(
                [
                    (0.5 * (bin_edges[i] + bin_edges[i + 1]))
                    * theHist[i]
                    * (bin_edges[i + 1] - bin_edges[i])
                    for i in range(len(theHist))
                ]
            )
            std68 = np.sqrt(
                np.sum(
                    [
                        ((0.5 * (bin_edges[i] + bin_edges[i + 1])) - mean) ** 2
                        * theHist[i]
                        * (bin_edges[i + 1] - bin_edges[i])
                        for i in range(len(theHist))
                    ]
                )
            )
            print("mean:", mean, "68% int.:", std68)
            return std68, mean - std68, mean + std68, MPV
        return 0.5 * (high - low), low, high, MPV

    def root_file_get_hist_and_edges(root_file, hist_name, rebin_factor=1):
        # Rebin_factor: if 1, no rebinning, if 2, combine every 2 bins, etc.
        # print("Available histograms:", [key.GetName() for key in root_file.GetListOfKeys()])
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
        if rebin_factor > 1:
            # rebin y and edges
            n_bins_rebinned = len(y) // rebin_factor
            y_rebinned = np.array(
                [
                    np.sum(y[i * rebin_factor : (i + 1) * rebin_factor])
                    for i in range(n_bins_rebinned)
                ]
            )
            edges_rebinned = np.array(
                [edges[i * rebin_factor] for i in range(n_bins_rebinned)]
                + [edges[n_bins_rebinned * rebin_factor]]
            )
            return y_rebinned, edges_rebinned
        return y, edges

    bin_mid_points = []
    lo_hi_MPV = []
    sigmaEoverE = []
    responses = []
    bins_to_histograms = {}

    def is_twojet_proc(procname):
        return "vvbb" in procname or "vvqq" in procname or "vvgg" in procname

    for i in range(len(bins) - 1):
        hist_name = f"{root_histogram_prefix}{suffix}{neg_format(bins[i])}_{neg_format(bins[i+1])}"
        field_names.append(hist_name)
        rf = 1
        # if bins[i] == 0 and bins[i+1] == 25 and is_twojet_proc(procname):
        #    rf = 2  # to reduce statistical fluctuations when stats are low
        y, edges = root_file_get_hist_and_edges(f, hist_name, rebin_factor=rf)
        if y is None:
            print(f"Skipping bin [{bins[i]}, {bins[i+1]}] due to missing histogram")
            continue
        # print("y:", y, "edges:", edges)
        # plot the current bin histogram using y and edges (NORMALIZED)
        # ax_hist.step(edges[:-1], y, where="post", label=f"[{bins[i]}, {bins[i+1]}] GeV")
        bin_widths = np.diff(edges)
        area = np.sum(y * bin_widths)
        n_jets_in_bin = np.sum(y)
        if area != 0:
            y_normalized = y / area
        else:
            y_normalized = y
        ax_hist[0].step(
            edges[:-1],
            y_normalized,
            where="post",
            label=f"[{bins[i]}, {bins[i + 1]}] GeV (N={int(n_jets_in_bin)})",
        )
        ax_hist[1].step(
            edges[:-1],
            y_normalized,
            where="post",
            label=f"[{bins[i]}, {bins[i + 1]}] GeV (N={int(n_jets_in_bin)})",
        )
        ax_hist[2].step(
            edges[:-1],
            y_normalized,
            where="post",
            label=f"[{bins[i]}, {bins[i + 1]}] GeV (N={int(n_jets_in_bin)})",
        )
        bins_to_histograms[i] = [y_normalized, edges]
        yc = copy(y)
        if sigma_method == "std68":
            std68, low, high, MPV = get_std68(
                y, edges, percentage=0.683, epsilon=eps, wmin=wmin, wmax=wmax
            )
        elif sigma_method == "RMS":
            MPV = 0.5 * (edges[np.argmax(y)] + edges[np.argmax(y) + 1])
            mean = np.sum(
                [
                    (0.5 * (edges[i] + edges[i + 1]))
                    * yc[i]
                    * (edges[i + 1] - edges[i])
                    for i in range(len(yc))
                ]
            ) / np.sum([yc[i] * (edges[i + 1] - edges[i]) for i in range(len(yc))])
            RMS = np.sqrt(
                np.sum(
                    [
                        ((0.5 * (edges[i] + edges[i + 1])) - mean) ** 2
                        * yc[i]
                        * (edges[i + 1] - edges[i])
                        for i in range(len(yc))
                    ]
                )
                / np.sum([yc[i] * (edges[i + 1] - edges[i]) for i in range(len(yc))])
            )
            std68 = RMS
            low = mean - std68
            high = mean + std68
        elif sigma_method == "interquantile_range":
            MPV = 0.5 * (edges[np.argmax(y)] + edges[np.argmax(y) + 1])
            # compute the cumulative distribution
            s = np.sum(yc * np.diff(edges))
            if s != 0:
                yc = yc / s  # normalize the histogram to 1
            cumulative = np.cumsum(yc * np.diff(edges))
            # find the 15.85% and 84.15% quantiles
            low_idx = np.searchsorted(cumulative, 0.1585)
            high_idx = np.searchsorted(cumulative, 0.8415)
            low = edges[low_idx]
            high = edges[high_idx]
            std68 = 0.5 * (high - low)
        elif sigma_method == "DSCB":
            centers = 0.5 * (edges[1:] + edges[:-1])
            MPV_guess = 0.5 * (edges[np.argmax(y)] + edges[np.argmax(y) + 1])
            sigma_guess = (
                np.std(np.repeat(centers, yc.astype(int))) if np.sum(yc) > 0 else 1.0
            )
            p0 = [MPV_guess, sigma_guess, 1.5, 3.0, 1.5, 3.0, max(yc)]
            try:
                popt, _ = curve_fit(
                    double_crystal_ball, centers, yc, p0=p0, maxfev=10000
                )
                mu, sigma, alphaL, nL, alphaR, nR, norm = popt
                MPV = mu
                std68 = sigma
                low, high = mu - sigma, mu + sigma
            except RuntimeError:
                print("⚠️ DSCB fit failed; reverting to RMS.")
                return get_result_for_process(y, edges, sigma_method="RMS")
        elif sigma_method == "gaussian_fit":
            # Fit a Gaussian to this histogram from 0.85 to 1.15
            centers = 0.5 * (edges[1:] + edges[:-1])
            mask = (centers >= 0.75) & (centers <= 1.15)
            centers_fit = centers[mask]
            yc_fit = yc[mask]
            if len(centers_fit) < 3:
                print("⚠️ Not enough points to fit a Gaussian; reverting to RMS.")
                return get_result_for_process(y, edges, sigma_method="RMS")

            def gaussian(x, mu, sigma, norm):
                return norm * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

            # fig, ax = plt.subplots()
            # make a tiny plot of what we are actually fitting and save as fit.pdf
            # ax.step(centers_fit, yc_fit, where="mid", label="Data")
            # fig.savefig("fit.pdf")
            mean_guess = 0.5 * (
                centers_fit[np.argmax(yc_fit)] + centers_fit[np.argmax(yc_fit) + 1]
            )
            sigma_guess = (
                np.std(np.repeat(centers_fit, yc_fit.astype(int)))
                if np.sum(yc_fit) > 0
                else 1.0
            )
            print("Mean guess:", mean_guess, "Sigma guess:", sigma_guess)
            p0 = [mean_guess, sigma_guess, max(yc_fit)]
            try:
                popt, _ = curve_fit(gaussian, centers_fit, yc_fit, p0=p0, maxfev=10000)
                mu, sigma, norm = popt
                MPV = mu
                std68 = sigma
                low, high = mu - sigma, mu + sigma
            except RuntimeError:
                print("⚠️ Gaussian fit failed; reverting to RMS.")
                return get_result_for_process(y, edges, sigma_method="RMS")
        else:
            raise ValueError(f"Unknown sigma method: {sigma_method}")
        lo_hi_MPV.append([low, high, MPV])
        bin_mid = 0.5 * (bins[i] + bins[i + 1])
        if (not np.isnan(bin_mid)) and (not np.isnan(std68)) and n_jets_in_bin > 50000:
            # More than 10k statistics for a good fit with the fine binning that we are using
            bin_mid_points.append(bin_mid)
            if divide_by_MPV:
                sigmaEoverE.append(std68 / MPV)
            else:
                sigmaEoverE.append(std68)
            responses.append(MPV)
            print(
                f"Bin [{bins[i]}, {bins[i+1]}]: {method} = {std68:.4f}, low = {low:.4f}, high = {high:.4f}, MPV={MPV}, N={np.sum(yc)} N_in_bin={n_jets_in_bin}"
            )
        else:
            print("NaN encountered in bin mid-point calculation / not enough statistics. Number of jets in bin:", n_jets_in_bin)
    ax_hist[0].legend(fontsize=9)
    ax_hist[0].set_xlabel(r"$E_{reco} / E_{true}$")
    ax_hist[0].set_ylabel("Entries")
    ax_hist[1].set_xlabel(r"$E_{reco} / E_{true}$")
    ax_hist[1].set_ylabel("Entries")
    ax_hist[1].set_xlim([wmin, wmax])
    ax_hist[2].set_xlim([wmin, wmax])

    ax_hist[2].set_yscale("log")
    ax_hist[2].set_xlabel(r"$E_{reco} / E_{true}$")
    ax_hist[2].set_ylabel("Entries")
    return (
        bin_mid_points,
        sigmaEoverE,
        fig_hist,
        responses,
        bins_to_histograms,
        lo_hi_MPV,
        field_names,
    )


bin_to_histograms_storage = {}
bin_to_histograms_storage_neutral = {}
method_low_high_mid_point_storage = {}


def get_func_fit(mid_points, Rs, confusion_term=True, min_E=0.0):
    # basically filter out points with mid_points < min_E
    mask = np.array(mid_points) >= min_E
    mid_points = np.array(mid_points)[mask]
    Rs = np.array(Rs)[mask]
    if confusion_term:

        def resolution_func(E, a, b, c):
            return np.sqrt((a / np.sqrt(E)) ** 2 + b**2 + (c / E) ** 2)

        popt, pcov = curve_fit(
            resolution_func, mid_points, Rs, p0=[0.5, 0.03, 0.1], maxfev=10000
        )
        # return a smooth function throughout mid_points
    else:

        def resolution_func(E, a, b):
            return np.sqrt((a / np.sqrt(E)) ** 2 + b**2)

        popt, pcov = curve_fit(
            resolution_func, mid_points, Rs, p0=[0.5, 0.03], maxfev=10000
        )
        # return a smooth function throughout mid_points
    xs = np.linspace(min(mid_points), max(mid_points), 100)
    ys = resolution_func(xs, *popt)
    popt = np.abs(popt)
    return xs, ys, popt, pcov


process_popt_storage = (
    {}
)  # Store the covariance matrix and optimal parameters for each process

method_to_color = {"std68": "blue", "gaussian_fit": "purple"}

jet_part_to_histogram_prefix = {
    "_charged": "binned_E_Charged_reco_over_true_FullGenJet_",
    "_neutral": "binned_E_Neutral_reco_over_true_FullGenJet_",
    "_all": "binned_E_reco_over_true_",
    "_photons": "binned_E_Photon_reco_over_true_FullGenJet_",
}

jet_parts_to_process = ["_all"] if args.angles_only else ["_photons", "_neutral", "_charged", "_all"]
for jet_part in jet_parts_to_process:
    if not args.angles_only:
        fig_resolution_per_process, ax_resolution_per_process = plt.subplots(
            len(processList), 2, figsize=(8, 4 * len(processList)), sharex=False
        )
        fig_resolution_per_process_Njets, ax_resolution_per_process_Njets = plt.subplots(
            2, 2, figsize=(9, 9), sharex=False
        )
    # Left column: resolution, right column: response. Comparison of gaussian_fit and std68
    for method in ["std68"]:
        print("-----------------------------------------------------------")
        print("Using peak width method:", method)
        # fig, ax = plt.subplots(2, 1, figsize=(8, 6))
        # the same as above but make it (10, 6) and make the upper plot 2/3 and lower plot 1/3 of the height
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
        for proc_idx, process in enumerate(sorted(list(processList.keys()))):
            if jet_part == "_all":
                (E_theta, sigma_theta, fig_theta_hist, response_theta, _, _, _) = get_result_for_process(
                    process,
                    sigma_method=method,
                    root_histogram_prefix="binned_deltaTheta_",
                    wmin=-0.05,
                    wmax=0.05,
                    divide_by_MPV=False
                )
                xs_theta, ys_theta, popt_theta, pcov_theta = get_func_fit(
                    E_theta, sigma_theta, confusion_term=False
                )
                (E_phi, sigma_phi, fig_phi_hist, response_phi, _, _, _) = get_result_for_process(
                    process,
                    sigma_method=method,
                    root_histogram_prefix="binned_deltaPhi_",
                    wmin=-0.05,
                    wmax=0.05,
                    divide_by_MPV=False
                )
                fig_theta_hist.savefig(os.path.join(outputDir, "bins_theta_{}.pdf".format(process)))
                fig_phi_hist.savefig(os.path.join(outputDir, "bins_phi_{}.pdf".format(process)))
                clr = PROCESS_COLORS.get(process, f"C{proc_idx}")
                xs_phi, ys_phi, popt_phi, pcov_phi = get_func_fit(
                    E_phi, sigma_phi, confusion_term=False
                )
                ax_theta[0].plot(E_theta, sigma_theta, "x", color=clr)
                ax_theta[0].plot(
                    xs_theta,
                    ys_theta,
                    LINE_STYLES[process],
                    color=clr,
                    label=HUMAN_READABLE_PROCESS_NAMES[process]
                    + f" {print_params(popt_theta)}",
                )
                ax_theta[1].plot(E_theta, response_theta, ".--", label=process, color=clr)
                ax_phi[0].plot(E_phi, sigma_phi, "x", color=clr)
                ax_phi[0].plot(
                    xs_phi,
                    ys_phi,
                    LINE_STYLES[process],
                    color=clr,
                    label=HUMAN_READABLE_PROCESS_NAMES[process]
                    + f" {print_params(popt_phi)}",
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
            if args.angles_only:
                continue
            (
                bin_mid_points,
                sigmaEoverE,
                fig_histograms,
                resp,
                bin_to_histograms,
                mpv_lo_hi,
                field_names,
            ) = get_result_for_process(
                process,
                sigma_method=method,
                root_histogram_prefix=jet_part_to_histogram_prefix[jet_part],
            )
            if process not in method_low_high_mid_point_storage:
                method_low_high_mid_point_storage[process] = {}
            if jet_part == "_all":
                method_low_high_mid_point_storage[process][method] = mpv_lo_hi
            if method == "std68" and jet_part == "_all":
                bin_to_histograms_storage[process] = bin_to_histograms
                fig_histograms.tight_layout()
                fig_histograms.savefig(
                    os.path.join(outputDir, "bins_{}_{}.pdf".format(process, method))
                )
            if method == "std68" and jet_part == "_neutral":
                bin_to_histograms_storage_neutral[process] = bin_to_histograms
                fig_histograms.tight_layout()
                fig_histograms.savefig(
                    os.path.join(
                        outputDir, "bins_NEUTRAL_{}_{}.pdf".format(process, method)
                    )
                )
            clr = PROCESS_COLORS.get(process, f"C{proc_idx}")
            if len(bin_mid_points) < 2:
                print(
                    f"Not enough points to fit for process {process} using method {method}. Skipping."
                )
                continue
            xs, ys, popt, pcov = get_func_fit(
                bin_mid_points, sigmaEoverE, confusion_term="CaloJet" in inputDir
            )
            print(f"Fitted parameters for {process} using {method}: {popt}")
            if process not in process_popt_storage:
                process_popt_storage[process] = {}
            process_popt_storage[process][method + jet_part] = (
                popt,
                pcov,
                xs,
                ys,
                bin_mid_points,
                sigmaEoverE,
                mpv_lo_hi,
                field_names,
            )
            ax[0].plot(bin_mid_points, sigmaEoverE, "x", color=clr)
            ax[0].plot(
                xs,
                ys,
                LINE_STYLES[process],
                color=clr,
                label=HUMAN_READABLE_PROCESS_NAMES[process] + f" {print_params(popt)}",
            )
            ax[1].plot(bin_mid_points, resp, ".--", label=process, color=clr)
            if method in method_to_color:
                ax_resolution_per_process[proc_idx, 0].plot(
                    bin_mid_points,
                    sigmaEoverE,
                    "x",
                    label=method + f" {print_params(popt)}",
                    color=method_to_color[method],
                )
                ax_resolution_per_process[proc_idx, 0].plot(
                    xs, ys, LINE_STYLES[process], color=method_to_color[method]
                )
                ax_resolution_per_process[proc_idx, 0].set_title(
                    HUMAN_READABLE_PROCESS_NAMES[process]
                )
                ax_resolution_per_process[proc_idx, 1].plot(
                    bin_mid_points,
                    resp,
                    ".--",
                    label=method,
                    color=method_to_color[method],
                )
            if method == "std68":
                if LINE_STYLES.get(process, "") == "-":  # full line
                    row = 0
                elif LINE_STYLES.get(process, "") == ":":  # dotted line
                    row = 1
                else:
                    row = None
                # Plot the b-jet containing events on the top plot, and light-jet containing events on the bottom plot
                if row is not None:
                    ax_resolution_per_process_Njets[row, 0].plot(
                        bin_mid_points,
                        sigmaEoverE,
                        "x",
                        label=HUMAN_READABLE_PROCESS_NAMES[process]
                        + f" {print_params(popt)}",
                        color=clr,
                    )
                    ax_resolution_per_process_Njets[row, 0].plot(
                        xs, ys, "--", color=clr
                    )
                    ax_resolution_per_process_Njets[row, 1].plot(
                        bin_mid_points,
                        resp,
                        ".--",
                        label=HUMAN_READABLE_PROCESS_NAMES[process],
                        color=clr,
                    )
            if not args.angles_only:
                ax_resolution_per_process[proc_idx, 0].set_xlabel("$E_{true}$ [GeV]")
                ax_resolution_per_process[proc_idx, 0].set_ylabel(r"$\sigma_E / E$")
                ax_resolution_per_process[proc_idx, 1].set_xlabel("$E_{true}$ [GeV]")
                ax_resolution_per_process[proc_idx, 1].set_ylabel("Response")
                # Also turn legend and grid on
                ax_resolution_per_process[proc_idx, 0].legend()
                ax_resolution_per_process[proc_idx, 0].grid(True)
                ax_resolution_per_process[proc_idx, 1].grid(True)
        if not args.angles_only:
            ax_resolution_per_process_Njets[0, 0].set_title("Final state containing b-jets")
            ax_resolution_per_process_Njets[1, 0].set_title(
                "Final state containing only light and gluon jets"
            )
            ax_resolution_per_process_Njets[0, 0].legend(
                title="l ∈ {u, d, s}; q ∈ {u, d, s, c, b}", fontsize=9.5, title_fontsize=8
            )
            ax_resolution_per_process_Njets[1, 0].legend(
                title="l ∈ {u, d, s}; q ∈ {u, d, s, c, b}", fontsize=9.5, title_fontsize=8
            )
            ax_resolution_per_process_Njets[0, 0].set_xlabel("$E_{true}$ [GeV]")
            ax_resolution_per_process_Njets[0, 1].set_xlabel("$E_{true}$ [GeV]")
            ax_resolution_per_process_Njets[0, 0].set_ylabel(r"$\sigma_E / E$")
            ax_resolution_per_process_Njets[1, 0].set_ylabel(r"$\sigma_E / E$")
            ax_resolution_per_process_Njets[1, 1].set_ylabel(r"$\sigma_E / E$")
            ax[0].legend()
            ax[0].set_xlabel("$E_{true}$ [GeV]")
            ax[0].set_ylabel(r"$\sigma_E / E$")
            ax[0].set_title(
                r"Jet Energy Resolution ($\frac{A}{\sqrt{E}}$ ⊕ \frac{B}{E} ⊕ $C$)"
            )
            ax[0].grid(True, alpha=0.3)
            ax[1].set_ylabel("Response")
            ax[1].set_xlabel("$E_{true}$ [GeV]")
            ax[1].grid()
            fig.tight_layout()
            fig.savefig(
                os.path.join(outputDir, "jet_energy_resolution_{}.pdf".format(method))
            )
        if jet_part == "_all":
            fig_theta.tight_layout()
            fig_theta.savefig(
                os.path.join(outputDir, "jet_theta_resolution_{}.pdf".format(method))
            )
            fig_phi.tight_layout()
            fig_phi.savefig(
                os.path.join(outputDir, "jet_phi_resolution_{}.pdf".format(method))
            )
    if not args.angles_only:
        fig_resolution_per_process.tight_layout()
        fig_resolution_per_process_Njets.tight_layout()
        fig_resolution_per_process.savefig(
            os.path.join(
                outputDir,
                "jet_energy_resolution_per_process_comparison{}.pdf".format(jet_part),
            )
        )
        fig_resolution_per_process_Njets.savefig(
            os.path.join(
                outputDir,
                "jet_energy_resolution_per_process_comparison_Njets{}.pdf".format(jet_part),
            )
        )

if not args.angles_only:
    pickle.dump(
        process_popt_storage,
        open(
            os.path.join(outputDir, "energy_fit_params_per_process.pkl"),
            "wb",
        ),
    )
    method_to_color = {
        "std68": "blue",
        "RMS": "orange",
        "interquantile_range": "green",
        "DSCB": "red",
        "gaussian_fit": "purple",
    }
else:
    import sys
    sys.exit(0)

### Plot each bin on a separate plot, but different processes on same plot ###
fig, ax = plt.subplots(
    len(binsE) - 1, 1, figsize=(6, 4 * (len(binsE) - 1)), sharex=False
)
fig_bins, ax_bins = plt.subplots(
    len(binsE) - 1, 1, figsize=(6, 4 * (len(binsE) - 1)), sharex=False
)

for i in range(len(binsE) - 1):
    for process in sorted(list(processList.keys())):
        y_normalized, edges = bin_to_histograms_storage[process][i]
        # plot on ax[i]
        bin_widths = np.diff(edges)
        ax[i].step(edges[:-1], y_normalized, where="post", label=process)
        ax_bins[i].step(edges[:-1], y_normalized, where="post", label=process)
        for method in method_low_high_mid_point_storage[process]:
            lo, hi, mpv = method_low_high_mid_point_storage[process][method][i]
            # plot vertical lines at lo, hi and mppv using method_to_color
            for _ax in [ax[i], ax_bins[i]]:
                _ax.axvline(
                    lo, color=method_to_color[method], linestyle="--", alpha=0.8
                )
                _ax.axvline(
                    hi, color=method_to_color[method], linestyle="--", alpha=0.8
                )
                _ax.axvline(
                    mpv, color=method_to_color[method], linestyle="-", alpha=0.8
                )
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

for method in ["std68"]:  # gaussian_fit optionally
    # fig, ax = plt.subplots(2, 1, figsize=(8, 6))
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]})
    for process in sorted(list(processList.keys())):
        bin_mid_points, sigmaEoverE, fig_histograms, resp, _, _, _ = (
            get_result_for_process(
                process, bins=bins_eta, suffix="eta_", sigma_method=method
            )
        )
        if method == "std68":
            fig_histograms.tight_layout()
            fig_histograms.savefig(
                os.path.join(outputDir, "bins_eta_{}.pdf".format(process))
            )
        ax[0].plot(
            bin_mid_points,
            sigmaEoverE,
            LINE_STYLES[process],
            label=HUMAN_READABLE_PROCESS_NAMES[process],
            color=PROCESS_COLORS[process],
        )
        ax[1].plot(
            bin_mid_points,
            resp,
            LINE_STYLES[process],
            label=HUMAN_READABLE_PROCESS_NAMES[process],
            color=PROCESS_COLORS[process],
        )
    ax[0].legend()
    ax[0].set_xlabel("$\eta$")
    ax[0].set_ylabel(r"$\sigma_E / E$")
    ax[0].set_title("Jet Energy Resolution vs $\Theta$")
    ax[0].grid(True, alpha=0.3)
    ax[1].grid(True, alpha=0.3)
    ax[1].set_title("Energy Response vs $\Theta$")
    ax[1].set_xlabel("$\eta$")
    ax[1].set_ylabel("Response")
    ax[1].set_xlabel("$\eta$")
    ax[1].grid()
    fig.tight_layout()
    fig.savefig(
        os.path.join(outputDir, "jet_Eta_resolution_data_points_{}.pdf".format(method))
    )

for method in ["std68"]:
    # fig, ax = plt.subplots(2, 1, figsize=(8, 6))
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={"height_ratios": [2, 1]})
    for process in sorted(list(processList.keys())):
        bin_mid_points, sigmaEoverE, fig_histograms, resp, _, _, _ = (
            get_result_for_process(
                process, bins=bins_costheta, suffix="costheta_", sigma_method=method
            )
        )
        if method == "std68":
            fig_histograms.tight_layout()
            fig_histograms.savefig(
                os.path.join(outputDir, "bins_CosTheta_{}.pdf".format(process))
            )
        ax[0].plot(bin_mid_points, sigmaEoverE, "x", color=PROCESS_COLORS[process])
        ax[1].plot(bin_mid_points, resp, "x", color=PROCESS_COLORS[process])
        ax[0].plot(
            bin_mid_points,
            sigmaEoverE,
            LINE_STYLES[process],
            label=HUMAN_READABLE_PROCESS_NAMES[process],
            color=PROCESS_COLORS[process],
        )
        ax[1].plot(
            bin_mid_points,
            resp,
            LINE_STYLES[process],
            label=HUMAN_READABLE_PROCESS_NAMES[process],
            color=PROCESS_COLORS[process],
        )
    ax[0].legend()
    ax[0].set_xlabel(r"cos $\theta$ [GeV]")
    ax[0].set_ylabel(r"$\sigma_E / E$")
    ax[0].set_title("Jet Energy Resolution vs Jet Energy")
    ax[0].grid(True, alpha=0.3)
    ax[1].grid(True, alpha=0.3)
    ax[1].set_title("Energy Resolution vs Jet Angle")
    ax[1].set_xlabel(r"cos $\theta$ [GeV]")
    ax[1].set_ylabel("Response")
    ax[1].set_xlabel(r"cos $\theta$ [GeV]")
    ax[1].grid()
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            outputDir, "jet_CosTheta_resolution_data_points_{}.pdf".format(method)
        )
    )
