"""RDataFrame graph for the per-event displays (src/event_displays.py).

Reproduces the selection chain that fills `h_mH_gen` - the "Physics" curve of
the mH decomposition figures - and then defines a flat, snapshot-able payload
per event: the stable gen particles, the gen jets, which jets the Higgs partons
matched to, and the invisible energy.

The point of the figure is to show the events that populate a given mH region,
so the chain here calls exactly the same helpers as `src/histmaker.py` rather
than re-deriving anything: `filter_MC_and_reco_particles`,
`compute_jets_from_args`, `get_hist_jet_eta_and_energy`,
`get_Higgs_mass_with_truth_matching`, and the same Filter strings.
"""
from src.histmaker_tools.jets import compute_jets_from_args
from src.histmaker_tools.jet_level_statistics import get_hist_jet_eta_and_energy
from src.histmaker_tools.particle_filter import filter_MC_and_reco_particles
from src.histmaker_tools.truth_matching import get_Higgs_mass_with_truth_matching

GEN_JET_VARIABLE = "GenJetFastJet"
RECO_JET_VARIABLE = "RecoJetFastJet"

# Columns AsNumpy() pulls out per selected event. Scalars first, then the jagged
# per-particle / per-jet / per-parton blocks.
SCALAR_COLUMNS = [
    "mH_gen",
    "mH_gen_all_jets",
    "mH_visible_truth",
    "mH_hard_partons",
    "E_nu_total",
    "E_nu_from_H",
    "E_vis_clustered",
    "n_gen_jets_total",
]
PARTICLE_COLUMNS = [
    "part_eta",
    "part_phi",
    "part_pt",
    "part_pdg",
    "part_energy",
    "part_px",
    "part_py",
    "part_pz",
    "part_jet_index",
]
JET_COLUMNS = [
    "jet_eta",
    "jet_phi",
    "jet_pt",
    "jet_energy",
    "jet_px",
    "jet_py",
    "jet_pz",
]
PARTON_COLUMNS = [
    "parton_eta",
    "parton_phi",
    "parton_pt",
    "parton_energy",
    "parton_pdg",
    "parton_to_jet",
]
PAYLOAD_COLUMNS = SCALAR_COLUMNS + PARTICLE_COLUMNS + JET_COLUMNS + PARTON_COLUMNS


def _serialize(df, prefix, collection, keys=("eta", "phi", "pt", "pdg", "energy")):
    """Unpack Utils::serialize_event(collection) into flat RVec columns.

    serialize_event returns (eta, phi, pt, pdg, mass, energy) as std::vectors;
    wrap each in an RVec so AsNumpy hands back something numpy can consume.
    """
    tuple_name = f"_ser_{prefix}"
    df = df.Define(tuple_name, f"FCCAnalyses::Utils::serialize_event({collection})")
    slots = {"eta": (0, "float"), "phi": (1, "float"), "pt": (2, "float"),
             "pdg": (3, "int"), "mass": (4, "float"), "energy": (5, "float")}
    for key in keys:
        index, ctype = slots[key]
        df = df.Define(
            f"{prefix}_{key}",
            f"ROOT::VecOps::RVec<{ctype}>(std::get<{index}>({tuple_name}))",
        )
    return df


def _momenta(df, prefix, collection):
    """px/py/pz of a Vec_rp, needed for the constituent-sum consistency check."""
    for axis in ("x", "y", "z"):
        df = df.Define(
            f"{prefix}_p{axis}",
            f"FCCAnalyses::ReconstructedParticle::get_p{axis}({collection})",
        )
    return df


def build_event_display_graph(df, args, n_jets, n_higgs_jets, apply_matched_filter=True):
    """Define every payload column. Returns (df, n_before_filters_node).

    args needs .jet_algorithm, .AK_radius, .energy_recovery, .ideal_matching and
    .jet_matching_radius - a types.SimpleNamespace is enough.
    """
    # --- the histmaker's selection chain, in the same order (histmaker.py:132-183)
    df = df.Define(
        "MC_part_idx",
        "FCCAnalyses::ZHfunctions::get_MC_quark_index_for_Higgs("
        "Particle, _Particle_daughters.index, false)",
    )
    df = df.Filter(f"MC_part_idx.size() == {n_higgs_jets}", "higgs partons found")
    df = filter_MC_and_reco_particles(df)
    df = compute_jets_from_args(df, args, n_jets)
    df, _ = get_hist_jet_eta_and_energy(df, RECO_JET_VARIABLE, GEN_JET_VARIABLE)
    df = df.Define(
        "reco_gen_jet_matching",
        "FCCAnalyses::JetTools::greedy_matching({}, {}, {})".format(
            RECO_JET_VARIABLE, GEN_JET_VARIABLE, args.jet_matching_radius
        ),
    )
    df = df.Define(
        "matched_genjet_E_and_all_genjet_E",
        "FCCAnalyses::JetTools::get_matched_and_all_E(reco_gen_jet_matching, {})".format(
            GEN_JET_VARIABLE
        ),
    )
    df = df.Define(
        "matched_genjet_energies", "get<0>(matched_genjet_E_and_all_genjet_E)"
    )
    if apply_matched_filter:
        # histmaker.py:179-183. On by default, because that is what the
        # PF_Durham tree behind the mH decomposition figures used.
        df = df.Filter(
            "(matched_genjet_energies.size() == genjet_energies.size()) && "
            "(genjet_energies.size() == {})".format(n_jets),
            "all jets matched",
        )
    df = get_Higgs_mass_with_truth_matching(
        df,
        genjets_field=GEN_JET_VARIABLE,
        recojets_field=RECO_JET_VARIABLE,
        expected_num_jets=n_higgs_jets,
        matching_radius=args.jet_matching_radius,
        n_jets=n_jets,
    )

    # --- the quantity the windows are cut on, plus context for the annotation box
    df = df.Define("mH_gen", "inv_mass_gen")
    df = df.Define("mH_gen_all_jets", "inv_mass_gen_all")
    df = df.Define("mH_visible_truth", "inv_mass_stable_gt_particles_from_higgs")
    df = df.Define("mH_hard_partons", "inv_mass_MC_part")
    df = df.Define("n_gen_jets_total", "(int) FastJet_jets.jets.size()")

    # --- what to draw
    df = _serialize(df, "part", "stable_gen_particles")
    df = _momenta(df, "part", "stable_gen_particles")
    df = _serialize(df, "jet", GEN_JET_VARIABLE, keys=("eta", "phi", "pt", "energy"))
    df = _momenta(df, "jet", GEN_JET_VARIABLE)
    df = _serialize(
        df, "parton", "MC_part_asjets", keys=("eta", "phi", "pt", "energy")
    )

    # MC_part_asjets goes through ZHfunctions::select_rp, which never sets PDG
    # (functions.h:606-632), so parton flavour has to come back from the MC record.
    df = df.Define(
        "parton_pdg",
        "FCCAnalyses::EventDisplay::parton_pdgs(Particle, MC_part_idx)",
    )
    # Index of the gen jet each Higgs parton matched to (-1 = unmatched).
    df = df.Define(
        "parton_to_jet", "ROOT::VecOps::RVec<int>(HardP_to_GenJet_mapping)"
    )
    df = df.Define(
        "part_jet_index",
        "FCCAnalyses::EventDisplay::constituent_jet_index("
        "FastJet_jets, (int) stable_gen_particles.size(), {})".format(n_jets),
    )

    # --- invisible / out-of-acceptance energy
    df = df.Define(
        "E_nu_total",
        "FCCAnalyses::EventDisplay::neutrino_energy(Particle, gt_labels, false)",
    )
    df = df.Define(
        "E_nu_from_H",
        "FCCAnalyses::EventDisplay::neutrino_energy(Particle, gt_labels, true)",
    )
    df = df.Define(
        "E_vis_clustered",
        "(float) ROOT::VecOps::Sum("
        "FCCAnalyses::ReconstructedParticle::get_e(stable_gen_particles))",
    )
    return df
