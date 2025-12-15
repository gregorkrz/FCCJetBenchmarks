def get_hist_jet_distances(df, RecoJetVariable, GenJetVariable):

    # Compute statistics for number of reconstructed jets, as well as distances between the jets

    df = df.Define("n_jets", "{}.size()".format(RecoJetVariable))
    df = df.Define("n_gen_jets", "{}.size()".format(GenJetVariable))
    df = df.Define(
        "distance_between_genjets",
        "FCCAnalyses::JetTools::get_pairwise_distances({})".format(GenJetVariable),
    )
    df = df.Define(
        "distance_between_recojets",
        "FCCAnalyses::JetTools::get_pairwise_distances({})".format(RecoJetVariable),
    )
    df = df.Define(
        "min_distance_between_genjets",
        "FCCAnalyses::Utils::min:with_default_value(FCCAnalyses::JetTools::get_jet_distances({}))".format(
            format(GenJetVariable)
        ),
    )
    df = df.Define(
        "min_distance_between_recojets",
        "FCCAnalyses::Utils::min:with_default_value(FCCAnalyses::ZHfunctions::get_jet_distances({}))".format(
            RecoJetVariable
        ),
    )
    # Will be different for each process with e+e- kt algorithm
    hist_njets = df.Histo1D(
        ("h_njets", "Number of reconstructed jets;N_jets;Events", 10, 0, 10), "n_jets"
    )
    hist_ngenjets = df.Histo1D(
        ("h_ngenjets", "Number of generated jets;N_genjets;Events", 10, 0, 10),
        "n_gen_jets",
    )
    hist_dist_jets_gen = df.Histo1D(
        (
            "h_dist_jets_gen",
            "Distance between gen jets;#DeltaR(jet_i, jet_j);Events",
            100,
            0,
            5,
        ),
        "distance_between_genjets",
    )
    hist_dist_jets_reco = df.Histo1D(
        (
            "h_dist_jets_reco",
            "Distance between reco jets;#DeltaR(jet_i, jet_j);Events",
            100,
            0,
            5,
        ),
        "distance_between_recojets",
    )
    hist_min_dist_jets_gen = df.Histo1D(
        (
            "h_min_dist_jets_gen",
            "Min distance between gen jets;min #DeltaR(jet_i, jet_j);Events",
            100,
            0,
            5,
        ),
        "min_distance_between_genjets",
    )
    hist_min_dist_jets_reco = df.Histo1D(
        (
            "h_min_dist_jets_reco",
            "Min distance between reco jets;min #DeltaR(jet_i, jet_j);Events",
            100,
            0,
            5,
        ),
        "min_distance_between_recojets",
    )
    histograms = [
        hist_njets,
        hist_ngenjets,
        hist_dist_jets_gen,
        hist_dist_jets_reco,
        hist_min_dist_jets_gen,
        hist_min_dist_jets_reco,
    ]
    return df, histograms


def get_hist_jet_eta_and_energy(df, RecoJetVariable, GenJetVariable):
    df = df.Define(
        "jet_etas", "FCCAnalyses::JetTools::get_jet_eta({})".format(RecoJetVariable)
    )
    df = df.Define(
        "genjet_etas",
        "FCCAnalyses::JetTools::get_jet_eta({})".format(GenJetVariable),
    )
    h_eta = df.Histo1D(("h_eta", "eta of reco jets;eta;Events", 100, -5, 5), "jet_etas")
    h_eta_gen = df.Histo1D(
        ("h_eta_gen", "eta of gen jets;eta;Events", 100, -5, 5), "genjet_etas"
    )
    df = df.Define(
        "jet_energies",
        "FCCAnalyses::JetTools::get_energy({})".format(RecoJetVariable),
    )
    df = df.Define(
        "genjet_energies",
        "FCCAnalyses::JetTools::get_energy({})".format(GenJetVariable),
    )
    h_E_jet = df.Histo1D(
        ("h_E_all_reco_jets", "E of reco jet;E_reco;Events", 100, 0, 300),
        "jet_energies",
    )
    h_E_genjet = df.Histo1D(
        ("h_E_all_gen_jets", "E of gen jet;E_gen;Events", 100, 0, 300),
        "genjet_energies",
    )
    return df, [h_eta, h_eta_gen, h_E_jet, h_E_genjet]
