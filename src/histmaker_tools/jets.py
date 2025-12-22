def _compute_anti_kt_jets(
    df, jet_clustering_function, energy_recovery, output_name, n_jets_energy_recovery
):
    """

    :param jet_clustering_function: Either JetClustering::clustering_antikt or JetClustering::clustering_ee_gen_antikt
    :return:
    """
    if energy_recovery:
        output_name_jets = "full_jets_{}".format(output_name)
    else:
        output_name_jets = output_name
    df = df.Define(output_name_jets, jet_clustering_function)
    if energy_recovery:
        df = df.Define(
            output_name,
            "FCCAnalyses::ZHfunctions::energy_recovery({}, {})".format(
                output_name_jets, n_jets_energy_recovery
            ),
        )
    return df


def get_jets(
    df,
    vec_rp_name,
    N_Durham=-1,
    AK_radius=-1,
    output_name="FastJet_jets",
    use_ee_AK=False,
    AK_energy_recovery=False,
    AK_energy_recovery_N_jets=2,
):
    """
    This function computes jets for a given collection of reconstructed particles or MC particles.
    vec_rp_name: name of the vector of ReconstructedParticles in the dataframe on which to perform jet clustering
    N_Durham: if set to > 0, jet clustering
    AK_radius: if set to > 0, anti-kt (use_ee_AK=False) or generalized e+e- anti-kt (use_ee_AK=True)
    name:
    """
    df = df.Define(
        "rp_px_{}".format(output_name),
        "FCCAnalyses::ReconstructedParticle::get_px({})".format(vec_rp_name),
    )
    df = df.Define(
        "rp_py_{}".format(output_name),
        "FCCAnalyses::ReconstructedParticle::get_py({})".format(vec_rp_name),
    )
    df = df.Define(
        "rp_pz_{}".format(output_name),
        "FCCAnalyses::ReconstructedParticle::get_pz({})".format(vec_rp_name),
    )
    df = df.Define(
        "rp_m_{}".format(output_name),
        "FCCAnalyses::ReconstructedParticle::get_mass({})".format(vec_rp_name),
    )
    df = df.Define(
        "fj_in_{}".format(output_name),
        "FCCAnalyses::JetClusteringUtils::set_pseudoJets_xyzm(rp_px_{},rp_py_{},rp_pz_{},rp_m_{})".format(
            output_name, output_name, output_name, output_name
        ),
    )
    if N_Durham > 0:
        print("Using Durham jet clustering algorithm with N={}".format(N_Durham))
        df = df.Define(
            output_name,
            "JetClustering::clustering_ee_kt(2, {}, 1, 0)(fj_in_{})".format(
                N_Durham, output_name
            ),
        )
    else:
        assert AK_radius > 0
        if not use_ee_AK:
            jets_func = (
                "JetClustering::clustering_antikt({}, 0, 0, 0, 0)(fj_in_{})".format(
                    AK_radius, output_name
                )
            )
            print("Using AK with R=", AK_radius)
        else:
            jets_func = (
                "JetClustering::clustering_ee_genkt({}, 0, 0, 0)(fj_in_{})".format(
                    AK_radius, output_name
                )
            )
            print("Using generalized e+e- AK with R=", AK_radius)
        df = _compute_anti_kt_jets(
            df,
            jets_func,
            AK_energy_recovery,
            output_name,
            n_jets_energy_recovery=AK_energy_recovery_N_jets,
        )
    return df


def get_jets_from_MC_reco_links(
    df, output_name="FastJet_jets", gen_jet_input="GenJet_jets"
):
    """
    Compute the jets from gen_jet_input
    """
    df = df.Define(
        output_name,
        "FCCAnalyses::ZHfunctions::match_genjet_constituents_to_reco_particles(ReconstructedParticlesEtaFilter, {}, mc2rp, stable_gen_particles_idx)".format(
            gen_jet_input
        ),
    )
    return df


def compute_jets_from_args(df, args, N_jets):
    # Compute the jets
    kwargs = {}
    if args.jet_algorithm == "Durham":
        kwargs["N_Durham"] = N_jets
    elif args.jet_algorithm in ["AK", "EEAK"]:
        kwargs["AK_radius"] = args.AK_radius
        assert args.AK_radius > 0
        if args.jet_algorithm == "EEAK":
            kwargs["use_ee_AK"] = True
        if args.energy_recovery:
            kwargs["AK_energy_recovery"] = True
            kwargs["AK_energy_recovery_N_jets"] = N_jets
    elif args.jet_algorithm == "CaloJetDurham":
        kwargs["N_Durham"] = N_jets
    else:
        raise ValueError("Unknown jet algorithm: {}".format(args.jet_algorithm))

    # Compute the Gen Jets
    df = get_jets(df, "stable_gen_particles", output_name="FastJet_jets", **kwargs)

    if args.ideal_matching:
        assert (
            args.jet_algorithm != "CaloJetDurham"
        ), "MC-reco links not implemented for CaloJets"
        df = get_jets_from_MC_reco_links(
            df, output_name="FastJet_jets_reco", gen_jet_input="FastJet_jets"
        )
    elif args.jet_algorithm == "CaloJetDurham":
        df = df.Define("FastJet_jets_reco", "CaloJetDurham")
    else:
        # Compute reco jets
        df = get_jets(
            df,
            "ReconstructedParticlesEtaFilter",
            output_name="FastJet_jets_reco",
            **kwargs,
        )
    first_k = N_jets
    df = df.Define(
        "GenJetFastJet",
        "FCCAnalyses::ZHfunctions::fastjet_to_vec_rp_jet(FastJet_jets, {})".format(
            first_k
        ),
    )
    if not args.jet_algorithm == "CaloJetDurham":
        df = df.Define(
            "RecoJetFastJet",
            "FCCAnalyses::ZHfunctions::fastjet_to_vec_rp_jet(FastJet_jets_reco, {})".format(
                first_k
            ),
        )
        # Store the neutral and charged components of the jets
        df = df.Define(
            "RecoJetFastJetNC",
            "FCCAnalyses::ZHfunctions::fastjet_to_vec_rp_jet_split_based_on_charge(FastJet_jets_reco, ReconstructedParticlesEtaFilter, {})".format(
                first_k
            ),
        )
        df = df.Define(
            "GenJetFastJetNC",
            "FCCAnalyses::ZHfunctions::fastjet_to_vec_rp_jet_split_based_on_charge(FastJet_jets, stable_gen_particles, {})".format(
                first_k
            ),
        )
    else:
        df = df.Define("RecoJetFastJet", "FastJet_jets_reco")
    return df
