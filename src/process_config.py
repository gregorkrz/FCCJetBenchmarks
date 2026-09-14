# This file contains information about the processes: The number of jets, the number of jets the Higgs decays into, as
# well as plotting metadata - colors, line styles, human-readable names of processes...
# It also holds the radius grid and the jet-algorithm family table, which used to be
# duplicated in generate_analysis_jobs.py, joint_plots.py and presentation_jer_plots.py.

# The radius grid shared by every radius-scan family. Single source of truth:
# scripts/generate_analysis_jobs.py, src/plotting/joint_plots.py and
# src/plotting/presentation_jer_plots.py all import it.
RADIUS_SCAN = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]


def radius_to_str(radius):
    """0.4 -> '04', 1.2 -> '12' - the suffix used in every method directory name."""
    return f"{int(round(radius * 10)):02d}"


# The radius-scan families, keyed by the name used for --algos and --family.
#
#   prefix         canonical method-directory prefix
#   legacy_prefix  the pre-2026-09 directory name, kept so an un-renamed tree still
#                  plots (and plots correctly labelled). See the note below.
#   label / suffix how the family is written in legends
#   exclusive_N    True  -> cluster to exactly NUMBER_OF_JETS[process] jets
#                  False -> inclusive clustering, jet multiplicity varies with R
#
# NAMING HISTORY: the directories called PF_AntiKtR* were never anti-kT. Until
# 2026-09, src/histmaker_tools/jets.py called
# JetClustering::clustering_ee_genkt with only four arguments, so the exponent
# fell back to the FCCAnalyses default of 0 - which is Cambridge/Aachen, not
# anti-kT (-1). Those histograms are valid e+e- C/A data; only the name was
# wrong, so they were relabelled and the directories renamed to
# PF_EECambridgeR*, with no re-run. See README, "Correction (2026-09)".
JET_FAMILIES = {
    "ee-ca": dict(
        prefix="PF_EECambridgeR",
        legacy_prefix="PF_AntiKtR",
        label="ee-C/A",
        suffix="",
        exclusive_N=False,
    ),
    "ee-ca-er": dict(
        prefix="PF_E_recovery_EECambridgeR",
        legacy_prefix="PF_E_recovery_AntiKtR",
        label="ee-C/A",
        suffix="-ER",
        exclusive_N=False,
    ),
    # e+e- kT (exponent +1), clustered inclusively so that R genuinely changes the
    # jets. Exclusive-to-N kT was measured to be bit-identical to Durham for
    # R >~ 0.6 (it is the R -> infinity limit of ee_genkt with p=+1), so an
    # exclusive scan over conventional radii would just reproduce PF_Durham.
    "ee-kt": dict(
        prefix="PF_EEKtR",
        legacy_prefix=None,
        label="ee-$k_T$",
        suffix="",
        exclusive_N=False,
    ),
    # Genuine e+e- anti-kT (exponent -1), first run 2026-09. NB the prefix is
    # PF_EEAntiKtR, which is deliberately distinct from the legacy PF_AntiKtR
    # (no "EE") used by the Cambridge/Aachen directories above.
    "ee-akt": dict(
        prefix="PF_EEAntiKtR",
        legacy_prefix=None,
        label="ee-anti-$k_T$",
        suffix="",
        exclusive_N=False,
    ),
}


def family_prefixes(family_key):
    """Canonical prefix first, then the legacy one if the family has one."""
    fam = JET_FAMILIES[family_key]
    return [p for p in (fam["prefix"], fam.get("legacy_prefix")) if p]


def method_radius(method_name, family_key=None):
    """Radius encoded in a method directory name, or None.

    'PF_EECambridgeR04' -> 0.4, 'PF_EEKtR12' -> 1.2, and the legacy
    'PF_AntiKtR04' -> 0.4. With family_key given, only that family's prefixes
    are tried, so PF_E_recovery_* is not mistaken for its plain twin.
    """
    keys = [family_key] if family_key else list(JET_FAMILIES)
    # Longest prefix first: PF_E_recovery_EECambridgeR must win over PF_EECambridgeR.
    prefixes = sorted(
        {p for k in keys for p in family_prefixes(k)}, key=len, reverse=True
    )
    for prefix in prefixes:
        if method_name.startswith(prefix):
            digits = method_name[len(prefix) :]
            if digits.isdigit():
                return int(digits) / 10
    return None


def method_family(method_name):
    """Family key a method directory belongs to, or None for Durham/CaloJets."""
    for key in sorted(
        JET_FAMILIES, key=lambda k: max(len(p) for p in family_prefixes(k)), reverse=True
    ):
        for prefix in family_prefixes(key):
            if method_name.startswith(prefix) and method_name[len(prefix) :].isdigit():
                return key
    return None


NUMBER_OF_JETS = {
    # 6 jets (blue)
    "p8_ee_ZH_6jet_ecm240": 6,
    "p8_ee_ZH_6jet_HF_ecm240": 6,
    "p8_ee_ZH_6jet_LF_ecm240": 6,
    # 4 jets (violet-magenta hues)
    "p8_ee_ZH_bbbb_ecm240": 4,
    "p8_ee_ZH_qqbb_ecm240": 4,
    "p8_ee_ZH_bbgg_ecm240": 4,
    "p8_ee_ZH_qqgg_ecm240": 4,
    "p8_ee_ZH_qqqq_ecm240": 4,
    # 2 jets (teal-green hues)
    "p8_ee_ZH_vvbb_ecm240": 2,
    "p8_ee_ZH_vvgg_ecm240": 2,
    "p8_ee_ZH_vvqq_ecm240": 2,
}

NUMBER_OF_HIGGS_JETS = {
    "p8_ee_ZH_qqbb_ecm240": 2,
    "p8_ee_ZH_6jet_ecm240": 4,
    "p8_ee_ZH_vvbb_ecm240": 2,
    "p8_ee_ZH_bbbb_ecm240": 2,
    "p8_ee_ZH_vvgg_ecm240": 2,
    "p8_ee_ZH_vvqq_ecm240": 2,
    "p8_ee_ZH_qqqq_ecm240": 2,
    "p8_ee_ZH_6jet_HF_ecm240": 4,
    "p8_ee_ZH_6jet_LF_ecm240": 4,
    "p8_ee_ZH_bbgg_ecm240": 2,
    "p8_ee_ZH_qqgg_ecm240": 2,
}

PROCESS_TO_ROW_COL = {
    # three columns for 2,4,6 jets
    "p8_ee_ZH_vvbb_ecm240": (0, 4),
    "p8_ee_ZH_vvgg_ecm240": (0, 1),
    "p8_ee_ZH_vvqq_ecm240": (0, 0),
    "p8_ee_ZH_bbbb_ecm240": (1, 4),
    "p8_ee_ZH_qqgg_ecm240": (1, 1),
    "p8_ee_ZH_qqbb_ecm240": (1, 2),
    "p8_ee_ZH_bbgg_ecm240": (1, 3),
    "p8_ee_ZH_qqqq_ecm240": (1, 0),
    "p8_ee_ZH_6jet_ecm240": (2, 1),
    "p8_ee_ZH_6jet_LF_ecm240": (2, 0),
    "p8_ee_ZH_6jet_HF_ecm240": (2, 2),
}

HUMAN_READABLE_PROCESS_NAMES = {
    "p8_ee_ZH_6jet_ecm240": "Z(→qq)H(→WW→qqqq) (all f.)",
    "p8_ee_ZH_6jet_HF_ecm240": "Z(→bb)H(→WW→bqbq)",
    "p8_ee_ZH_6jet_LF_ecm240": "Z(→qq)H(→WW→qqqq)",
    "p8_ee_ZH_bbbb_ecm240": "Z(→bb)H(→bb)",
    "p8_ee_ZH_qqbb_ecm240": "Z(→qq)H(→bb)",
    "p8_ee_ZH_bbgg_ecm240": "Z(→bb)H(→gg)",
    "p8_ee_ZH_qqgg_ecm240": "Z(→qq)H(→gg)",
    "p8_ee_ZH_qqqq_ecm240": "Z(→qq)H(→qq)",
    "p8_ee_ZH_vvbb_ecm240": "Z(→νν)H(→bb)",
    "p8_ee_ZH_vvgg_ecm240": "Z(→νν)H(→gg)",
    "p8_ee_ZH_vvqq_ecm240": "Z(→νν)H(→qq)",
}

LINE_STYLES = {
    # Containing only light-flavour jets: dotted; containing only b-jets: full line,
    # Containing mixture of some sorts or just gluons: dashed
    "p8_ee_ZH_6jet_ecm240": "--",
    "p8_ee_ZH_6jet_HF_ecm240": "-",
    "p8_ee_ZH_6jet_LF_ecm240": ":",
    "p8_ee_ZH_bbbb_ecm240": "-",
    "p8_ee_ZH_qqbb_ecm240": "--",
    "p8_ee_ZH_bbgg_ecm240": "--",
    "p8_ee_ZH_qqgg_ecm240": "--",
    "p8_ee_ZH_qqqq_ecm240": ":",
    "p8_ee_ZH_vvbb_ecm240": "-",
    "p8_ee_ZH_vvgg_ecm240": "--",
    "p8_ee_ZH_vvqq_ecm240": ":",
}

PROCESS_COLORS = {
    # 6 jets (blue)
    "p8_ee_ZH_6jet_ecm240": "#0067A5",
    "p8_ee_ZH_6jet_HF_ecm240": "#0082C8",
    "p8_ee_ZH_6jet_LF_ecm240": "#339EDD",
    # 4 jets (violet-magenta hues)
    "p8_ee_ZH_bbbb_ecm240": "#B832A0",
    "p8_ee_ZH_qqbb_ecm240": "#7A3CBF",
    "p8_ee_ZH_bbgg_ecm240": "#D45ECF",
    "p8_ee_ZH_qqgg_ecm240": "#D890E0",
    "p8_ee_ZH_qqqq_ecm240": "#E3B1F0",
    # 2 jets (teal-green hues)
    "p8_ee_ZH_vvbb_ecm240": "#1B9E77",
    "p8_ee_ZH_vvgg_ecm240": "#33AF8A",
    "p8_ee_ZH_vvqq_ecm240": "#7CCBA2",
}
