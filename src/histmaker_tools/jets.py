# FastJet's generalized e+e- algorithm is selected by its exponent p, the 6th
# argument of JetClustering::clustering_ee_genkt:
#
#     p = -1   e+e- anti-kT
#     p =  0   e+e- Cambridge/Aachen
#     p = +1   e+e- kT
#
# with d_ij = min(E_i^2p, E_j^2p) * (1 - cos(theta_ij)) / (1 - cos(R)) and
# d_iB = E_i^2p. The 1/(1 - cos R) factor is a global constant, so R never
# reorders the d_ij among themselves: the whole R dependence enters through the
# d_ij vs d_iB (beam) comparison. Consequently p = +1 run *exclusively* to N
# jets is bit-identical to Durham for R >~ 0.6 (Durham is the R -> infinity
# limit), which is why the ee-kT family below is clustered inclusively.
#
# HISTORY / BUG: until 2026-09 this file called clustering_ee_genkt with only
# FOUR arguments, so p silently fell back to the FCCAnalyses default of 0.
# Every histogram in a PF_AntiKtR* or PF_E_recovery_AntiKtR* directory is
# therefore e+e- Cambridge/Aachen, NOT anti-kT. Those runs are valid C/A data
# and were relabelled rather than re-run; the directories are now
# PF_EECambridgeR*. Always pass the exponent explicitly.
EE_GENKT_EXPONENT = {
    "EEAKT": -1.0,  # e+e- anti-kT (genuine; first run 2026-09)
    "EECA": 0.0,  # e+e- Cambridge/Aachen (what the historical EEAK runs did)
    "EEKT": 1.0,  # e+e- kT
}

# Accepted --jet-algorithm spellings that map onto clustering_ee_genkt, resolved
# to the canonical key of EE_GENKT_EXPONENT.
#
# CAREFUL: EEAK and EEAKT are NOT the same thing.
#   EEAK  - the historical name, kept because ~400 job scripts under jobs*/ use
#           it. It ran with no exponent argument, i.e. p = 0, so it resolves to
#           EECA (Cambridge/Aachen) and its output is unchanged.
#   EEAKT - genuine e+e- anti-kT, p = -1.
EE_GENKT_ALGORITHMS = {
    "EEAK": "EECA",
    "EECA": "EECA",
    "EEKT": "EEKT",
    "EEAKT": "EEAKT",
}

# Inclusive (False) or exclusive-to-N (True) clustering per algorithm. All three
# families are inclusive, so that R genuinely changes the jets: an exclusive-to-N
# kT scan would be degenerate, since exclusive ee_genkt with p=+1 reproduces
# Durham bit-for-bit for R >~ 0.6. (FastJet also warns that exclusive jets are
# only well defined for p >= 0, which rules it out for anti-kT anyway.)
EE_GENKT_EXCLUSIVE_N = {"EEAKT": False, "EECA": False, "EEKT": False}


def _define_jets_with_optional_energy_recovery(
    df, jet_clustering_function, energy_recovery, output_name, n_jets_energy_recovery
):
    """Define the jet collection, optionally merging surplus jets into the N leading ones.

    :param jet_clustering_function: a JetClustering::clustering_* Define string
    :param energy_recovery: if set, cluster into full_jets_<output_name> and run
        ZHfunctions::energy_recovery to end up with n_jets_energy_recovery jets
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
    ee_genkt_exponent=None,
    ee_genkt_exclusive_N=None,
    AK_energy_recovery=False,
    AK_energy_recovery_N_jets=2,
):
    """
    This function computes jets for a given collection of reconstructed particles or MC particles.

    vec_rp_name: name of the vector of ReconstructedParticles in the dataframe on which to perform jet clustering
    N_Durham: if > 0, exclusive ee-kT (Durham) clustering to exactly this many jets
    AK_radius: if > 0, radius-based clustering. Which algorithm is chosen by
        ee_genkt_exponent:
          None -> JetClustering::clustering_antikt, the hadron-collider anti-kT
          float -> JetClustering::clustering_ee_genkt with that exponent p
                   (-1 anti-kT, 0 Cambridge/Aachen, +1 kT; see EE_GENKT_EXPONENT)
    ee_genkt_exclusive_N: None -> inclusive clustering (multiplicity varies with R);
        an int -> cluster to exactly that many jets
    output_name: name of the FCCAnalysesJet column to define
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
        if ee_genkt_exponent is None:
            jets_func = (
                "JetClustering::clustering_antikt({}, 0, 0, 0, 0)(fj_in_{})".format(
                    AK_radius, output_name
                )
            )
            print("Using hadron-collider anti-kt with R=", AK_radius)
        else:
            if ee_genkt_exclusive_N is None:
                exclusive, cut = 0, 0.0  # inclusive, no pT cut
            else:
                exclusive, cut = 2, float(ee_genkt_exclusive_N)  # exactly N jets
            # (radius, exclusive, cut, sorted=0 pT-ordered, recombination=0 E-scheme, exponent)
            jets_func = (
                "JetClustering::clustering_ee_genkt({R}, {excl}, {cut}, 0, 0, {p})"
                "(fj_in_{name})".format(
                    R=AK_radius,
                    excl=exclusive,
                    cut=cut,
                    p=float(ee_genkt_exponent),
                    name=output_name,
                )
            )
            print(
                "Using generalized e+e- kt: p={} R={} exclusive={} cut={}".format(
                    float(ee_genkt_exponent), AK_radius, exclusive, cut
                )
            )
        df = _define_jets_with_optional_energy_recovery(
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
    elif args.jet_algorithm == "AK" or args.jet_algorithm in EE_GENKT_ALGORITHMS:
        kwargs["AK_radius"] = args.AK_radius
        assert args.AK_radius > 0
        if args.jet_algorithm in EE_GENKT_ALGORITHMS:
            # EEAK is the historical spelling and resolves to EECA: those runs
            # passed no exponent, so they were Cambridge/Aachen all along.
            canonical = EE_GENKT_ALGORITHMS[args.jet_algorithm]
            kwargs["ee_genkt_exponent"] = EE_GENKT_EXPONENT[canonical]
            if EE_GENKT_EXCLUSIVE_N[canonical]:
                kwargs["ee_genkt_exclusive_N"] = N_jets
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
    # In a plain reco run, also build the ideal-matching jets (gen-jet grouping
    # with the PFlow partners of each constituent). Having both jet definitions
    # in the *same* run is what makes it possible to require the Higgs jets to be
    # found by both, so the mH definitions can be compared on identical events.
    df = df.Define("_have_alt_jets", "true")
    if not args.ideal_matching and args.jet_algorithm != "CaloJetDurham":
        df = get_jets_from_MC_reco_links(
            df, output_name="FastJet_jets_ideal", gen_jet_input="FastJet_jets"
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
    if not args.ideal_matching and args.jet_algorithm != "CaloJetDurham":
        df = df.Define(
            "IdealJetFastJet",
            "FCCAnalyses::ZHfunctions::fastjet_to_vec_rp_jet(FastJet_jets_ideal, {})".format(
                first_k
            ),
        )
    return df
