def get_Higgs_mass_with_truth_matching(
    df,
    genjets_field="GenJetDurhamN4",
    recojets_field="JetDurhamN4",
    expected_num_jets=-1,
    matching_radius=0.3,
):
    """
    Get mH by matching the genjets to
    :param df:
    :param genjets_field:
    :param recojets_field:
    :param expected_num_jets:
    :param matching_radius:
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
    df = df.Define(
        "inv_mass_gen_all",
        "FCCAnalyses::ZHfunctions::invariant_mass({})".format(genjets_field),
    )
    df = df.Define(
        "inv_mass_reco_all",
        "FCCAnalyses::ZHfunctions::invariant_mass({})".format(recojets_field),
    )
    return df
