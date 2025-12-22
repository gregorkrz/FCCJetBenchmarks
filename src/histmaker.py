# Arguments: input_dir, output_dir, jet_algo, ak_radius (for AK), keep_only_fully_matched_events, ideal_matching
import os
import sys
import pickle
import argparse

from src.histmaker_tools.truth_matching import get_Higgs_mass_with_truth_matching
from src.process_config import NUMBER_OF_JETS, NUMBER_OF_HIGGS_JETS
from src.histmaker_tools.jets import compute_jets_from_args
from src.histmaker_tools.particle_filter import filter_MC_and_reco_particles
from src.histmaker_tools.jet_level_statistics import get_hist_jet_distances
from src.histmaker_tools.hit_level_statistics import get_calo_hit_statistics
from src.histmaker_tools.binning import bin_quantity
from src.histmaker_tools.jet_level_statistics import get_hist_jet_eta_and_energy

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="Input directory containing the dataset")
parser.add_argument(
    "--output", type=str, help="Output directory to store the produced files"
)
parser.add_argument(
    "--jet-algorithm",
    type=str,
    choices=["Durham", "AK", "EEAK", "CaloJetDurham"],
    default="Durham",
    help="Jet algorithm to use",
)
parser.add_argument(
    "--AK-radius", type=float, default=0.6, help="The radius parameter for AK and ee-AK"
)
parser.add_argument(
    "--jet-matching-radius",
    type=float,
    default=1.0,
    help="The radius parameter for AK and ee-AK",
)
parser.add_argument(
    "--no-filter-fully-matched",
    action="store_true",
    default=False,
    help="By default, we only keep the events where all the jets are fully matched. If toggled, this filter will be turned off.",
)
parser.add_argument(
    "--ideal-matching",
    action="store_true",
    default=False,
    help="If toggled, the reco jets will be computed using the MC-reco links from the gen jets.",
)
parser.add_argument(
    "--only-dataset",
    type=str,
    default="",
    help="If set, it will only process one dataset.",
)

parser.add_argument(
    "--energy-recovery",
    action="store_true",
    help="If set, energy recovery will be used for the anti-kt or generalized e+e- anti-kt jet algorithm",
)

if "--" in sys.argv:
    argv_after_sep = sys.argv[sys.argv.index("--") + 1 :]
    args = parser.parse_args(argv_after_sep)
else:
    args = parser.parse_args()

inputDir = args.input
print("Using input dir:", inputDir)

frac = 1
processList = {}

for folder in os.listdir(inputDir):
    procname = folder.replace(".root", "")
    if args.only_dataset != "" and args.only_dataset != procname:
        continue
    processList[procname] = {"fraction": frac}
    if procname not in NUMBER_OF_HIGGS_JETS or procname not in NUMBER_OF_JETS:
        raise Exception(f"Process {procname} not known")
print("ProcList:", processList)
bins_jetE = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
bins_eta = [-5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 5]
bins_costheta = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]

procDict = "FCCee_procDict_winter2023_IDEA.json"

includePaths = [
    "histmaker_functions/jet_tools.h",
    "histmaker_functions/utils.h",
    "histmaker_functions/functions.h",
]

outputDir = args.output
os.makedirs(outputDir, exist_ok=True)

Fully_Matched_Only = not args.no_filter_fully_matched

if Fully_Matched_Only:
    print(
        "Keeping only fully matched events in the output histograms! Jet matching efficiency will be 1."
    )


GenJetVariable = "GenJetFastJet"
RecoJetVariable = "RecoJetFastJet"
nCPUS = -1
doScale = False
intLumi = 5000000  # 5 /ab

bins_count_jets = (5, 0, 5)


def point_format(number):
    return str(number).replace(".", "p")


def neg_format(number):
    # put n5 for -5
    if number < 0:
        return point_format("n{}".format(abs(number)))
    else:
        return point_format(number)


# build_graph function that contains the analysis logic, cuts and histograms (mandatory)
def build_graph(df, dataset):
    print("############## Doing dataset:", dataset, "##############")
    histograms = []  # output histogram root files
    df = df.Define(
        "MC_part_idx",
        "FCCAnalyses::ZHfunctions::get_MC_quark_index_for_Higgs(Particle, _Particle_daughters.index, false)",
    )
    # Sometimes, the W+W- are not decayed in the ZH -> 6 jet events. We filter out these events (<1%)
    df = df.Filter("MC_part_idx.size() == {}".format(NUMBER_OF_HIGGS_JETS[dataset]))

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    df, calo_hit_histograms = get_calo_hit_statistics(df)
    histograms += calo_hit_histograms

    # Create Gen Particle and Reco Particle collections
    df = filter_MC_and_reco_particles(df)

    # Compute the jets from arguments
    df = compute_jets_from_args(df, args, NUMBER_OF_JETS[dataset])
    df, histograms_jet_eta_E = get_hist_jet_eta_and_energy(
        df, RecoJetVariable, GenJetVariable
    )
    histograms += histograms_jet_eta_E

    df = df.Define(
        "reco_gen_jet_matching",
        "FCCAnalyses::JetTools::greedy_matching({}, {}, {})".format(
            RecoJetVariable, GenJetVariable, args.jet_matching_radius
        ),
    )
    df, hist_jet_distances = get_hist_jet_distances(df, RecoJetVariable, GenJetVariable)
    histograms += hist_jet_distances
    df = df.Define(
        "matched_genjet_E_and_all_genjet_E",
        "FCCAnalyses::JetTools::get_matched_and_all_E(reco_gen_jet_matching, {})".format(
            GenJetVariable
        ),
    )
    df = df.Define(
        "matched_genjet_energies", "get<0>(matched_genjet_E_and_all_genjet_E)"
    )
    if Fully_Matched_Only:
        # Filter by matched_genjet_energies.size() == genjet_energies.size() &&
        # genjet_energies.size() == nJets_processList[dataset]
        print(
            "Filtering! Number of matched jets should be equal to",
            NUMBER_OF_JETS[dataset],
        )
        df = df.Filter(
            "(matched_genjet_energies.size() == genjet_energies.size()) && (genjet_energies.size() == {})".format(
                NUMBER_OF_JETS[dataset]
            )
        )
    weightsum_after_filtering = df.Sum("weight").GetValue()
    # Save to pickle: before and after filtering weightsum
    out_file = os.path.join(outputDir, f"basic_stats_{dataset}.pkl")
    with open(out_file, "wb") as f:  # Keep track of how many events were filtered out
        pickle.dump(
            {
                "before_filtering": weightsum.GetValue(),
                "after_filtering": weightsum_after_filtering,
            },
            f,
        )
    df = df.Define("all_genjet_energies", "get<1>(matched_genjet_E_and_all_genjet_E)")
    hist_genjet_all_energies = df.Histo1D(
        ("h_genjet_all_energies", "E of all gen jets;E_gen;Events", 10, 0, 200),
        "all_genjet_energies",
    )
    hist_genjet_matched_energies = df.Histo1D(
        ("h_genjet_matched_energies", "E of matched gen jets;E_gen;Events", 10, 0, 200),
        "matched_genjet_energies",
    )
    df = df.Define(
        "matching_processing",
        "FCCAnalyses::JetTools::get_E_reco_over_E_true(reco_gen_jet_matching, {}, {})".format(
            RecoJetVariable, GenJetVariable
        ),
    )
    if not args.jet_algorithm == "CaloJetDurham":
        # Use the RecoJetVariable+NC and GenJetVariable+NC for getting the charged and neutral components of the jets
        # First element of the pairs: neutral (N), second: charged (C)
        df = df.Define(
            "matching_processing_Neutral_part",
            "FCCAnalyses::ZHfunctions::get_energy_ratios_for_matched_jets(reco_gen_jet_matching, get<0>(RecoJetFastJetNC), get<0>(GenJetFastJetNC))",
        )
        df = df.Define(
            "jet_E_reco_over_E_true_Neutral_part",
            "std::get<0>(matching_processing_Neutral_part)",
        )
        df = df.Define(
            "matching_processing_Charged_part",
            "FCCAnalyses::ZHfunctions::get_energy_ratios_for_matched_jets(reco_gen_jet_matching, get<1>(RecoJetFastJetNC), get<1>(GenJetFastJetNC))",
        )
        df = df.Define(
            "jet_E_reco_over_E_true_Charged_part",
            "std::get<0>(matching_processing_Charged_part)",
        )
        df = df.Define(
            "matching_processing_Photon_part",
            "FCCAnalyses::ZHfunctions::get_energy_ratios_for_matched_jets(reco_gen_jet_matching, get<2>(RecoJetFastJetNC), get<2>(GenJetFastJetNC))",
        )
        df = df.Define(
            "jet_E_reco_over_E_true_Photon_part",
            "std::get<0>(matching_processing_Photon_part)",
        )
        df = df.Define(
            "genjet_Neutral_energies_matched",
            "std::get<2>(matching_processing_Neutral_part)",
        )
        df = df.Define(
            "genjet_Charged_energies_matched",
            "std::get<2>(matching_processing_Charged_part)",
        )
        df = df.Define(
            "genjet_Photon_energies_matched",
            "std::get<2>(matching_processing_Photon_part)",
        )
        df = df.Define(
            "matching_processing_Neutral_part_wrt_Full_GenJets",
            "FCCAnalyses::ZHfunctions::get_energy_ratios_for_matched_jets(reco_gen_jet_matching, get<0>(RecoJetFastJetNC), {})".format(
                GenJetVariable
            ),
        )
        df = df.Define(
            "jet_E_reco_over_E_true_Neutral_part_wrt_Full_GenJets",
            "std::get<0>(matching_processing_Neutral_part_wrt_Full_GenJets)",
        )
        df = df.Define(
            "matching_processing_Charged_part_wrt_Full_GenJets",
            "FCCAnalyses::ZHfunctions::get_energy_ratios_for_matched_jets(reco_gen_jet_matching, get<1>(RecoJetFastJetNC), {})".format(
                GenJetVariable
            ),
        )
        df = df.Define(
            "jet_E_reco_over_E_true_Charged_part_wrt_Full_GenJets",
            "std::get<0>(matching_processing_Charged_part_wrt_Full_GenJets)",
        )
        df = df.Define(
            "matching_processing_Photon_part_wrt_Full_GenJets",
            "FCCAnalyses::ZHfunctions::get_energy_ratios_for_matched_jets(reco_gen_jet_matching, get<2>(RecoJetFastJetNC), {})".format(
                GenJetVariable
            ),
        )
        df = df.Define(
            "jet_E_reco_over_E_true_Photon_part_wrt_Full_GenJets",
            "std::get<0>(matching_processing_Photon_part_wrt_Full_GenJets)",
        )
        df = df.Define(
            "genjet_Neutral_energies_matched_wrt_Full_GenJets",
            "std::get<2>(matching_processing_Neutral_part_wrt_Full_GenJets)",
        )
        df = df.Define(
            "genjet_Charged_energies_matched_wrt_Full_GenJets",
            "std::get<2>(matching_processing_Charged_part_wrt_Full_GenJets)",
        )
        df = df.Define(
            "genjet_Photon_energies_matched_wrt_Full_GenJets",
            "std::get<2>(matching_processing_Photon_part_wrt_Full_GenJets)",
        )
        # Print some of jet_E_reco_over_E_true_Neutral_part, jet_E_reco_over_E_true_Neutral_part
    df = df.Define("jet_E_reco_over_E_true", "std::get<0>(matching_processing)")
    # Print the first 5 elements
    df = df.Define("E_of_unmatched_reco_jets", "std::get<1>(matching_processing)")
    df = df.Define("num_unmatched_reco_jets", "E_of_unmatched_reco_jets.size()")
    df = df.Define("genjet_energies_matched", "std::get<2>(matching_processing)")
    df = df.Define("genjet_etas_matched", "std::get<3>(matching_processing)")
    df = df.Define(
        "genjet_costhetas_matched",
        "FCCAnalyses::Utils::get_costheta_from_eta(genjet_etas_matched)",
    )
    df = df.Define("num_matched_reco_jets", "genjet_energies_matched.size()")
    # Bin the jet_E_reco_over_E_true according to genjet_energies (bins [0, 50, 100, 150, 200])
    histograms += [
        hist_genjet_all_energies,
        hist_genjet_matched_energies,
    ]
    # this can be replaced with df, hist_E = bin_quantity(df, "
    df, hist_jetE = bin_quantity(
        df,
        "jet_E_reco_over_E_true",
        "genjet_energies_matched",
        bins=bins_jetE,
        output_prefix="binned_E_reco_over_true",
    )
    histograms += hist_jetE
    if args.jet_algorithm != "CaloJetDurham":
        # For studying the NH, charged, and photon components separately
        df, hist_jetE_neutral = bin_quantity(
            df,
            "jet_E_reco_over_E_true_Neutral_part",
            "genjet_energies_matched",
            bins=bins_jetE,
            output_prefix="binned_E_Neutral_reco_over_true_FullGenJet",
        )
        histograms += hist_jetE_neutral
        df, hist_jetE_charged = bin_quantity(
            df,
            "jet_E_reco_over_E_true_Charged_part",
            "genjet_energies_matched",
            bins=bins_jetE,
            output_prefix="binned_E_Charged_reco_over_true_FullGenJet",
        )
        histograms += hist_jetE_charged
        df, hist_jetE_photon = bin_quantity(
            df,
            "jet_E_reco_over_E_true_Photon_part",
            "genjet_energies_matched",
            bins=bins_jetE,
            output_prefix="binned_E_Photon_reco_over_true_FullGenJet",
        )
        histograms += hist_jetE_photon

    # For each energy and eta bin, count the number of unmatched reco jets over the number of gen jets in that bin,
    # and save these variables in the output file too
    df, histograms_E_binned_by_eta = bin_quantity(
        df,
        "jet_E_reco_over_E_true",
        "genjet_etas_matched",
        bins=bins_eta,
        output_prefix="binned_E_reco_over_true_eta",
    )
    histograms += histograms_E_binned_by_eta
    df, histograms_E_binned_by_cos_Theta = bin_quantity(
        df,
        "jet_E_reco_over_E_true",
        "genjet_costhetas_matched",
        bins=bins_costheta,
        output_prefix="binned_E_reco_over_true_costheta",
    )
    histograms += histograms_E_binned_by_cos_Theta
    h_jet_E_reco_over_E_true = df.Histo1D(
        (
            "h_jet_E_reco_over_E_true",
            "E_reco/E_true (fancy matching);E_reco / E_true;Events",
            300,
            0.5,
            1.5,
        ),
        "jet_E_reco_over_E_true",
    )
    # df = df.Define(
    #    "jet_E_reco_over_E_true_higheta",
    #    "FCCAnalyses::ZHfunctions::cut_by_quantity(jet_E_reco_over_E_true, genjet_etas_matched, -999999, -0.9)",
    # )
    histograms += [h_jet_E_reco_over_E_true]
    # Make a histogram of jet energies

    h_unmatched_reco_jets = df.Histo1D(
        (
            "h_unmatched_reco_jets",
            "E of unmatched reco jets;E_reco;Events",
            100,
            0,
            300,
        ),
        "E_of_unmatched_reco_jets",
    )
    histograms += [h_unmatched_reco_jets]

    df = get_Higgs_mass_with_truth_matching(
        df,
        genjets_field=GenJetVariable,
        recojets_field=RecoJetVariable,
        expected_num_jets=NUMBER_OF_HIGGS_JETS[dataset],
        matching_radius=args.jet_matching_radius,
    )
    df = df.Define(
        "matching_reco_with_partons",
        "FCCAnalyses::JetTools::greedy_matching({}, {}, {})".format(
            RecoJetVariable, "MC_part_asjets", args.jet_matching_radius
        ),
    )
    df = df.Define(
        "matching_proc_with_partons",
        "FCCAnalyses::ZHfunctions::get_energy_ratios_for_matched_jets(matching_reco_with_partons, {}, {})".format(
            RecoJetVariable, "MC_part_asjets"
        ),
    )
    df = df.Define(
        "ratio_jet_energies_matching_with_partons",
        "get<0>(matching_proc_with_partons)",
    )
    df = df.Define(
        "E_of_unmatched_reco_jets_with_partons",
        "get<1>(matching_proc_with_partons)",
    )
    h_ratio_matching_with_partons = df.Histo1D(
        (
            "h_ratio_matching_with_partons",
            "E_reco/E_parton;E_reco / E_parton;Events",
            300,
            0.5,
            1.5,
        ),
        "ratio_jet_energies_matching_with_partons",
    )
    df = df.Define(
        "inv_mass_all_gen_particles",
        "FCCAnalyses::ZHfunctions::invariant_mass(stable_gen_particles);",
    )

    h_mH_all_stable_part = df.Histo1D(
        (
            "h_mH_all_stable_part",
            "Invariant mass of all particles; Minv; Events",
            1000,
            0,
            250,
        ),
        "inv_mass_all_gen_particles",
    )
    # do the invariant mass of all reco particles
    df = df.Define(
        "inv_mass_all_reco_particles",
        "FCCAnalyses::ZHfunctions::invariant_mass(ReconstructedParticlesEtaFilter);",
    )
    hist_m_all_reco_particles = df.Histo1D(
        (
            "hist_calo_hist_E",
            "Invariant mass of all reco particles; Minv; Events",
            1000,
            0,
            250,
        ),
        "inv_mass_all_reco_particles",
    )
    h_mH_reco = df.Histo1D(
        ("h_mH_reco", "Higgs mass from reco jets;M_H (reco jets);Events", 1000, 0, 250),
        "inv_mass_reco",
    )
    h_mH_gen = df.Histo1D(
        ("h_mH_gen", "Higgs mass from gen jets;M_H (gen jets);Events", 1000, 0, 250),
        "inv_mass_gen",
    )
    h_mH_gen_all = df.Histo1D(
        (
            "h_mH_gen_all",
            "Higgs mass from all gen jets;M_H (all gen jets);Events",
            1000,
            0,
            250,
        ),
        "inv_mass_gen_all",
    )
    h_mH_reco_all = df.Histo1D(
        (
            "h_mH_reco_all",
            "Higgs mass from all reco jets;M_H (all reco jets);Events",
            1000,
            0,
            250,
        ),
        "inv_mass_reco_all",
    )
    histograms += [
        h_mH_reco,
        h_mH_gen,
        h_mH_gen_all,
        h_mH_reco_all,
        h_mH_all_stable_part,
        hist_m_all_reco_particles,
    ]
    h_mH_stable_gt_particles = df.Histo1D(
        (
            "h_mH_stable_gt_particles",
            "Higgs mass from stable gt particles;M_H (stable gt particles);Events",
            1000,
            0,
            250,
        ),
        "inv_mass_stable_gt_particles_from_higgs",
    )
    h_mH_reco_particles_matched = df.Histo1D(
        (
            "h_mH_reco_particles_matched",
            "Higgs mass from reco particles matched;M_H (reco particles matched);Events",
            1000,
            0,
            250,
        ),
        "inv_mass_reco_particles_matched_from_higgs",
    )
    h_mH_MC_part = df.Histo1D(
        (
            "h_mH_MC_part",
            "Higgs mass from initial MC part.;M_H (MC part);Events",
            1000,
            0,
            250,
        ),
        "inv_mass_MC_part",
    )
    histograms += [
        h_mH_stable_gt_particles,
        h_mH_reco_particles_matched,
        h_mH_MC_part,
        h_ratio_matching_with_partons,
    ]
    # Define a constant that is the df length after filtering and also save that
    return histograms, weightsum
