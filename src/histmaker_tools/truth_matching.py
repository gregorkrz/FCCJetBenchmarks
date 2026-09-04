def get_Higgs_mass_with_truth_matching(
    df,
    genjets_field="GenJetDurhamN4",
    recojets_field="JetDurhamN4",
    expected_num_jets=-1,
    matching_radius=0.3,
    recojets_fastjet_field=None,
    n_jets=-1,
    alt_recojets_field=None,
):
    """
    Get mH by matching the genjets to
    :param df:
    :param genjets_field:
    :param recojets_field:
    :param expected_num_jets:
    :param matching_radius:
    :param recojets_fastjet_field: the JetClustering::FCCAnalysesJet the reco
        jets came from (constituents needed for the provenance-based Higgs-jet
        selection). None -> that variable is not defined (e.g. CaloJets, whose
        jets are read from a branch and have no constituents).
    :param n_jets: number of jets kept per event (first_k in the pT sort)
    :param alt_recojets_field: a second reco-jet collection (the ideal-matching
        jets) to compute mH from in the same run, so that events can be required
        to have their Higgs jets found by *both* definitions
    :return:
    """
    df = df.Define(
        "MC_quark_idx", "FCCAnalyses::ZHfunctions::get_MC_quark_index(Particle)"
    )
    df = df.Define(
        "gt_labels",
        "FCCAnalyses::ZHfunctions::getGTLabels(MC_part_idx, Particle, _Particle_daughters.index);",
    )
    df = df.Define(
        "_gt_particles_from_higgs",
        "FCCAnalyses::ZHfunctions::select_gt_particles(gt_labels, Particle)",
    )
    df = df.Define("gt_part_from_H_idx", "_gt_particles_from_higgs.first")
    df = df.Define("stable_gt_particles_from_higgs", "_gt_particles_from_higgs.second")
    df = df.Define(
        "inv_mass_stable_gt_particles_from_higgs",
        "FCCAnalyses::ZHfunctions::invariant_mass(stable_gt_particles_from_higgs)",
    )
    df = df.Define(
        "reco_particles_matched_from_higgs",
        "FCCAnalyses::ZHfunctions::get_particles_from_mc2rp(gt_part_from_H_idx, mc2rp, ReconstructedParticlesEtaFilter)",
    )
    df = df.Define(
        "inv_mass_reco_particles_matched_from_higgs",
        "FCCAnalyses::ZHfunctions::invariant_mass(reco_particles_matched_from_higgs)",
    )
    # MC_part_asjets contains
    df = df.Define(
        "MC_part_asjets",
        "FCCAnalyses::ZHfunctions::select_rp(FCCAnalyses::ZHfunctions::vec_mc_to_rp(Particle), MC_part_idx)",
    )

    df = df.Define(
        "inv_mass_MC_part", "FCCAnalyses::ZHfunctions::invariant_mass(MC_part_asjets)"
    )
    df = df.Define(
        "HardP_to_GenJet_mapping",
        "FCCAnalyses::JetTools::greedy_matching(MC_part_asjets, {}, {})".format(
            genjets_field, matching_radius
        ),
    )
    df = df.Define(
        "HardP_to_RecoJet_mapping",
        "FCCAnalyses::ZHfunctions::merge_mappings(HardP_to_GenJet_mapping, reco_gen_jet_matching)",
    )
    df = df.Define(
        "filtered_jets",
        "FCCAnalyses::ZHfunctions::filter_jets({}, HardP_to_RecoJet_mapping)".format(
            recojets_field
        ),
    )
    # The mapping above composes parton->genjet with reco_gen_jet_matching, which
    # is a reco->gen map, so it only lands on the right reco jet where that
    # permutation is its own inverse (identity or pure pairwise swaps). A cycle of
    # length >= 3 picks a wrong (but same-size) jet subset, so the event survives
    # with a wrong mass. Redo it with a genuine gen->reco matching; greedy
    # matching isn't symmetric, so run it in that direction rather than inverting.
    df = df.Define(
        "gen_reco_jet_matching",
        "FCCAnalyses::JetTools::greedy_matching({}, {}, {})".format(
            genjets_field, recojets_field, matching_radius
        ),
    )
    df = df.Define(
        "HardP_to_RecoJet_mapping_fixed",
        "FCCAnalyses::ZHfunctions::merge_mappings(HardP_to_GenJet_mapping, gen_reco_jet_matching)",
    )
    df = df.Define(
        "filtered_jets_fixed",
        "FCCAnalyses::ZHfunctions::filter_jets({}, HardP_to_RecoJet_mapping_fixed)".format(
            recojets_field
        ),
    )
    df = df.Define(
        "inv_mass_reco_fixed",
        "FCCAnalyses::ZHfunctions::invariant_mass(filtered_jets_fixed, {})".format(
            expected_num_jets
        ),
    )
    df = df.Define(
        "filtered_jets_gen",
        "FCCAnalyses::ZHfunctions::filter_jets({}, HardP_to_GenJet_mapping)".format(
            genjets_field
        ),
    )
    df = df.Define(
        "inv_mass_reco",
        "FCCAnalyses::ZHfunctions::invariant_mass(filtered_jets, {})".format(
            expected_num_jets
        ),
    )
    df = df.Define(
        "inv_mass_gen",
        "FCCAnalyses::ZHfunctions::invariant_mass(filtered_jets_gen, {})".format(
            expected_num_jets
        ),
    )
    # Higgs jets chosen by particle provenance rather than by an angular match:
    # each jet is ranked by the energy fraction of its constituents whose MC
    # partner descends from a Higgs parton, and the top expected_num_jets are
    # summed. This always yields the expected number of jets, so it separates
    # particle-level misassignment from gen-reco jet-matching failures.
    if recojets_fastjet_field is not None:
        df = df.Define(
            "higgs_jets_by_provenance",
            "FCCAnalyses::ZHfunctions::select_higgs_jets_by_provenance("
            "{}, ReconstructedParticlesEtaFilter, rp2mc, gt_labels, {}, {})".format(
                recojets_fastjet_field, n_jets, expected_num_jets
            ),
        )
        df = df.Define(
            "inv_mass_reco_provenance",
            "FCCAnalyses::ZHfunctions::invariant_mass(higgs_jets_by_provenance, {})".format(
                expected_num_jets
            ),
        )
    # Same Higgs-jet selection applied to a second jet collection. Uses its own
    # gen->reco matching, so the two definitions are treated identically.
    if alt_recojets_field is not None:
        df = df.Define(
            "gen_altreco_jet_matching",
            "FCCAnalyses::JetTools::greedy_matching({}, {}, {})".format(
                genjets_field, alt_recojets_field, matching_radius
            ),
        )
        df = df.Define(
            "HardP_to_AltRecoJet_mapping",
            "FCCAnalyses::ZHfunctions::merge_mappings(HardP_to_GenJet_mapping, gen_altreco_jet_matching)",
        )
        df = df.Define(
            "filtered_jets_alt",
            "FCCAnalyses::ZHfunctions::filter_jets({}, HardP_to_AltRecoJet_mapping)".format(
                alt_recojets_field
            ),
        )
        df = df.Define(
            "inv_mass_reco_alt",
            "FCCAnalyses::ZHfunctions::invariant_mass(filtered_jets_alt, {})".format(
                expected_num_jets
            ),
        )
    df = df.Define(
        "inv_mass_gen_all",
        "FCCAnalyses::ZHfunctions::invariant_mass({})".format(genjets_field),
    )
    df = df.Define(
        "inv_mass_reco_all",
        "FCCAnalyses::ZHfunctions::invariant_mass({})".format(recojets_field),
    )
    return df
