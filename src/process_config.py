# This file contains information about the processes: The number of jets, the number of jets the Higgs decays into, as
# well as plotting metadata - colors, line styles, human-readable names of processes...

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
    ## For the ee-AK comparison
    "Calo Durham": "Calo Durham",
    "PF Durham": "PF Durham",
    "AK4": "ee-AK4",
    "AK6": "ee-AK6",
    "AK8": "ee-AK8",
    "AK10": "ee-AK10",
    "AK12": "ee-AK12",
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
    "Calo Durham": ":",
    "PF Durham": "-",
    "AK4": "--",
    "AK6": "--",
    "AK8": "--",
    "AK10": "--",
    "AK12": "--",
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
    # For the comparison of jet algorithms
    "Calo Durham": "#1f77b4",
    "PF Durham": "#1fa39c",
    "AK4": "#b99cf2",
    "AK6": "#a67ee8",
    "AK8": "#8c5cd9",
    "AK10": "#6e3ebf",
    "AK12": "#542ba2",
}
