import ROOT
import os
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", type=str, required=True)

if "--" in sys.argv:
    argv_after_sep = sys.argv[sys.argv.index("--") + 1 :]
    args = parser.parse_args(argv_after_sep)
else:
    args = parser.parse_args()

intLumi = 1.0
intLumiLabel = ""
# ana_tex = os.environ["FOLDER_NAME"]
ana_tex = ""

delphesVersion = "3.4.2"
energy = 240.0
collider = "FCC-ee"
formats = ["png", "pdf"]

inputDir = args.inputDir
outdir = os.path.join(inputDir, "plots_debug")

os.mkdir(outdir) if not os.path.exists(outdir) else None

plotStatUnc = True

colors = {}
color_presets = [
    ROOT.kRed,
    ROOT.kBlue + 1,
    ROOT.kCyan + 2,
    ROOT.kMagenta,
    ROOT.kOrange + 7,
    ROOT.kGray + 1,
    ROOT.kGreen,
    ROOT.kOrange,
    ROOT.kViolet,
    ROOT.kAzure + 1,
    ROOT.kPink + 1,
    ROOT.kTeal + 1,
]

procs = {"signal": {}, "backgrounds": {}}

legend = {}

i = 0

for file in sorted(os.listdir(inputDir)):
    if file.endswith(".root"):
        proc_name = file.replace(".root", "")
        procs["signal"][proc_name] = [proc_name]
        legend[proc_name] = proc_name
        colors[proc_name] = color_presets[i]
        i += 1

print("Procs:", procs, "Legend:", legend)

hists = {}

hists["h_jet_E_reco_over_E_true"] = {
    "output": "jet_E_reco_over_E_true",
    "logy": False,
    "stack": False,
    "ymin": 0,
    "xtitle": "E_reco/E_true (deltaR matching)",
    "ytitle": "Events",
}

hists["h_unmatched_reco_jets"] = {
    "output": "E_of_unmatched_reco_jets",
    "logy": False,
    "stack": False,
    # "rebin":    100,
    # "xmin":     120,
    # "xmax":     140,
    ##"ymin":     0,
    "xtitle": "E of unmatched reco jets",
    "ytitle": "Events",
}

hists["h_unmatched_reco_jets"] = {
    "output": "E_of_unmatched_reco_jets",
    "logy": False,
    "stack": False,
    # "rebin":    100,
    # "xmin":     120,
    # "xmax":     140,
    ##"ymin":     0,
    "xtitle": "E of unmatched reco jets",
    "ytitle": "Events",
}

hists["h_eta"] = {
    "output": "h_eta",
    "logy": False,
    "stack": False,
    # "rebin":    100,
    # "xmin":     120,
    # "xmax":     140,
    ##"ymin":     0,
    "xtitle": "eta of reco jets",
    "ytitle": "Events",
}

hists["h_eta_gen"] = {
    "output": "h_eta_gen",
    "logy": False,
    "stack": False,
    # "rebin":    100,
    # "xmin":     120,
    # "xmax":     140,
    ##"ymin":     0,
    "xtitle": "eta of gen jets",
    "ytitle": "Events",
}


hists["h_dist_jets_gen"] = {
    "output": "h_dist_jets_gen",
    "logy": False,
    "stack": False,
    "xtitle": "DeltaR between gen jets",
    "ytitle": "Pairs of gen jets",
}


hists["h_min_dist_jets_gen"] = {
    "output": "h_min_dist_jets_gen",
    "logy": False,
    "stack": False,
    "xtitle": "Min. DeltaR between gen jets",
    "ytitle": "Events",
}

hists["h_min_dist_jets_reco"] = {
    "output": "h_min_dist_jets_reco",
    "logy": False,
    "stack": False,
    "xtitle": "Min. DeltaR between reco jets",
    "ytitle": "Events",
}

hists["h_mH_reco"] = {
    "output": "h_mH_reco",
    "logy": False,
    "stack": False,
    "ymax": 20000,
    "xtitle": "Reconstructed Higgs mass (all matched jets)",
    "ytitle": "Events",
    # Normalize to 1
}

hists["h_mH_gen"] = {
    "output": "h_mH_gen",
    "logy": False,
    "stack": False,
    "ymax": 120000,
    "xtitle": " Higgs mass from gen particles (all matched jets)",
    "ytitle": "Events",
}


hists["h_njets"] = {
    "output": "h_njets",
    "logy": False,
    "stack": False,
    # "rebin":    100,
    # "xmin":     120,
    # "xmax":     140,
    ##"ymin":     0,
    "xtitle": "Number of reco jets",
    "ytitle": "Events",
}

hists["h_ngenjets"] = {
    "output": "h_ngenjets",
    "logy": False,
    "stack": False,
    # "rebin":    100,
    # "xmin":     120,
    # "xmax":     140,
    ##"ymin":     0,
    "xtitle": "Number of gen jets",
    "ytitle": "Events",
}

hists["h_mH_gen_all"] = {
    "output": "h_mH_gen_all",
    "logy": False,
    "stack": False,
    "xtitle": " Higgs mass from gen particles (all gen jets)",
    "ytitle": "Events",
}

hists["h_mH_reco_all"] = {
    "output": "h_mH_reco_all",
    "logy": False,
    "stack": False,
    "xtitle": "Reconstructed Higgs mass (all reco jets)",
    "ytitle": "Events",
}

hists["h_mH_all_stable_part"] = {
    "output": "h_mH_all_stable_part",
    "logy": False,
    "stack": False,
    "xtitle": "invariant mass of all particles",
    "ytitle": "Events",
}

hists["h_E_all_reco_jets"] = {
    "output": "h_E_all_reco_jets",
    "logy": False,
    "stack": False,
    "xtitle": "reco jet E",
    "ytitle": "Events",
}

hists["h_E_all_gen_jets"] = {
    "output": "h_E_all_genjets",
    "logy": True,
    "stack": False,
    "xtitle": "gen jet E",
    "ytitle": "Events",
}


hists["h_mH_stable_gt_particles"] = {
    "output": "h_mH_stable_gt_particles",
    "logy": True,
    "stack": False,
    "xtitle": "Higgs mass from stable GT particles",
    "ytitle": "Events",
}

hists["h_mH_reco_particles_matched"] = {
    "output": "h_mH_reco_particles_matched",
    "logy": True,
    "stack": False,
    "xtitle": "Higgs mass from reco particles matched from Higgs",
    "ytitle": "Events",
}

hists["h_mH_MC_part"] = {
    "output": "h_mH_MC_part",
    "logy": True,
    "stack": False,
    "xtitle": "Higgs mass from MC partons",
    "ytitle": "Events",
}
