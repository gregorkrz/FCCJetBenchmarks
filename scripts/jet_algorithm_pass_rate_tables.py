#!/usr/bin/env python3
# Print common-sample pass-rate tables for Durham vs inclusive e+e- algos.
# Usage: python3 scripts/jet_algorithm_pass_rate_tables.py <HIST_DIR> [--r 08] [--per-n]
import argparse
import os
import sys

import ROOT  # noqa: E402
ROOT.gROOT.SetBatch(True)

PROCS = [
    ("p8_ee_ZH_vvbb_ecm240", 2, r"$Z(\to\nu\bar\nu)\,H(\to b\bar b)$",         "Z(νν)H(bb)"),
    ("p8_ee_ZH_vvgg_ecm240", 2, r"$Z(\to\nu\bar\nu)\,H(\to gg)$",              "Z(νν)H(gg)"),
    ("p8_ee_ZH_vvqq_ecm240", 2, r"$Z(\to\nu\bar\nu)\,H(\to q\bar q)$",         "Z(νν)H(qq)"),
    ("p8_ee_ZH_bbbb_ecm240", 4, r"$Z(\to b\bar b)\,H(\to b\bar b)$",           "Z(bb)H(bb)"),
    ("p8_ee_ZH_bbgg_ecm240", 4, r"$Z(\to b\bar b)\,H(\to gg)$",                "Z(bb)H(gg)"),
    ("p8_ee_ZH_qqbb_ecm240", 4, r"$Z(\to q\bar q)\,H(\to b\bar b)$",           "Z(qq)H(bb)"),
    ("p8_ee_ZH_qqgg_ecm240", 4, r"$Z(\to q\bar q)\,H(\to gg)$",                "Z(qq)H(gg)"),
    ("p8_ee_ZH_qqqq_ecm240", 4, r"$Z(\to q\bar q)\,H(\to q\bar q)$",           "Z(qq)H(qq)"),
    ("p8_ee_ZH_6jet_HF_ecm240", 6, r"$Z(\to b\bar b)\,H(\to WW^*\to bqbq)$",         "Z(bb)H(WW→bqbq)"),
    ("p8_ee_ZH_6jet_LF_ecm240", 6, r"$Z(\to q\bar q)\,H(\to WW^*\to qqqq)$ (light f.)", "Z(qq)H(WW) LF"),
    ("p8_ee_ZH_6jet_ecm240",    6, r"$Z(\to q\bar q)\,H(\to WW^*\to qqqq)$ (all f.)",   "Z(qq)H(WW) all"),
]

ORDER = ["Durham", "antikT", "CA", "kT", "antikT_ER", "CA_ER", "kT_ER"]
HEADER = ["Durham", "anti-kT", "C/A", "kT", "anti-kT +ER", "C/A +ER", "kT +ER"]


def methods_for(r):
    return {
        "Durham":     "PF_Durham",
        "antikT":     f"PF_EEAntiKtR{r}",
        "CA":         f"PF_EECambridgeR{r}",
        "kT":         f"PF_EEKtR{r}",
        "antikT_ER":  f"PF_E_recovery_EEAntiKtR{r}",
        "CA_ER":      f"PF_E_recovery_EECambridgeR{r}",
        "kT_ER":      f"PF_E_recovery_EEKtR{r}",
    }


STAGE_HISTS = {
    "physics":  "h_mH_gen",          # gen-jet clustering finds all Higgs partons
    "det_phys": "h_mH_reco",         # reco jets find them via gen->reco mapping
    "common":   "h_mH_common_reco",  # both reco and ideal matching succeed
}


def _count_pass(f, hname):
    h = f.Get(hname)
    if not h:
        return 0
    return int(round(h.Integral(2, h.GetNbinsX() + 1)))


def pass_pcts(hist_dir, method_dir, proc, min_tot=1000):
    path = os.path.join(hist_dir, method_dir, f"{proc}.root")
    if not os.path.exists(path):
        return {k: None for k in STAGE_HISTS}
    f = ROOT.TFile.Open(path)
    tot_h = f.Get("h_mH_stable_gt_particles")
    tot = int(tot_h.GetEntries()) if tot_h else 0
    if tot < min_tot:
        f.Close()
        return {k: None for k in STAGE_HISTS}
    out = {k: 100.0 * _count_pass(f, h) / tot for k, h in STAGE_HISTS.items()}
    f.Close()
    return out


def pass_pct(hist_dir, method_dir, proc, min_tot=1000):
    return pass_pcts(hist_dir, method_dir, proc, min_tot)["common"]


def fmt(v, best, bold=True):
    if v is None:
        return "---"
    s = f"{v:.1f}"
    if bold and abs(v - best) < 0.05:
        return f"**{s}**"
    return s


def collect(hist_dir, r_for_proc):
    rows = []
    for proc, nj, _, short in PROCS:
        r = r_for_proc(nj)
        m = methods_for(r)
        vals = {k: pass_pct(hist_dir, m[k], proc) for k in ORDER}
        rows.append((nj, short, r, vals))
    return rows


def md_table(rows, title, note):
    out = [f"### {title}", "", note, ""]
    out.append("| Process | jets | R | " + " | ".join(HEADER) + " |")
    out.append("| --- | --- | --- | " + " | ".join(["---"] * len(HEADER)) + " |")
    last_n = None
    for nj, short, r, vals in rows:
        if last_n is not None and nj != last_n:
            out.append("|   |   |   | " + " | ".join([" "] * len(HEADER)) + " |")
        vv = [v for v in vals.values() if v is not None]
        best = max(vv) if vv else 0
        cells = [fmt(vals[k], best) for k in ORDER]
        out.append(f"| {short} | {nj} | {int(r)/10:.1f} | " + " | ".join(cells) + " |")
        last_n = nj
    out.append("")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hist_dir")
    ap.add_argument("--r", default="08", help="R for the fixed-R table (default 08)")
    ap.add_argument("--per-n", action="store_true",
                    help="Also print the per-N table (2j/4j at R=1.4, 6j at R=1.0)")
    ap.add_argument("--stages", action="store_true",
                    help="Also print Physics / Det+Phys / Common pass rates per method")
    args = ap.parse_args()

    r = args.r
    r_str = f"R={int(r)/10:.1f}"
    print(f"# Jet-algorithm pass-rate comparison (R={int(r)/10:.1f})")
    print()
    print(f"Dataset: `{args.hist_dir}`")
    print()
    print("Pass rate = 100 × events with mH_common_reco > 0 / all generated events "
          "(common sample: both standard and alt mapping succeed). Best in bold.")
    print()
    rows_fixed = collect(args.hist_dir, lambda nj: r)
    print(md_table(rows_fixed, f"Common sample, {r_str} for inclusive algos",
                   "Durham is exclusive-N; inclusive algos use the R above."))

    if args.per_n:
        rows_pern = collect(args.hist_dir, lambda nj: {2: "14", 4: "14", 6: "10"}[nj])
        print()
        print(md_table(rows_pern,
                       "Common sample, per-N R (2j & 4j: R=1.4, 6j: R=1.0)",
                       "6-jet rows at R=1.4 crash (can't form 6 jets); R=1.0 used instead."))

    if args.stages:
        print()
        print(f"## Pass rate by stage at R={int(r)/10:.1f}: Physics / Det+Phys / Common sample")
        print()
        print("Physics = gen-jet clustering alone finds all Higgs partons. "
              "Det+Phys = reco jets find them via gen→reco mapping. "
              "Common = both reco and ideal matching succeed (the mH-decomp sample). "
              "Best in each row in bold.")
        STAGE_LABEL = [("physics", "Physics"), ("det_phys", "Det+Phys"), ("common", "Common")]
        for stage_key, stage_name in STAGE_LABEL:
            print()
            print(f"### {stage_name}")
            print()
            print("| Process | jets | " + " | ".join(HEADER) + " |")
            print("| --- | --- | " + " | ".join(["---"] * len(HEADER)) + " |")
            last_n = None
            m = methods_for(r)
            for proc, nj, _, short in PROCS:
                if last_n is not None and nj != last_n:
                    print("|   |   | " + " | ".join([" "] * len(HEADER)) + " |")
                vals = {k: pass_pcts(args.hist_dir, m[k], proc)[stage_key] for k in ORDER}
                vv = [v for v in vals.values() if v is not None]
                best = max(vv) if vv else 0
                cells = [fmt(vals[k], best) for k in ORDER]
                print(f"| {short} | {nj} | " + " | ".join(cells) + " |")
                last_n = nj


if __name__ == "__main__":
    main()
