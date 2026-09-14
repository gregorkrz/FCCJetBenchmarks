"""Stage 1 of the event displays: pick events in three mH windows and dump them.

The mH decomposition figures show a broad "Physics" curve - `h_mH_gen`, the
Higgs mass from *gen* jets built out of stable gen particles, i.e. a perfect
detector, so its width is entirely a jet-definition / jet-assignment effect.
This script selects events far below, at, and far above that curve's peak and
writes out everything needed to draw them, so the tails can be inspected event
by event.

Nothing under $PATH_TO_HISTOGRAMS can be used for this: those files hold only
TH1Ds. The clustering has to be redone from the source EDM4hep files, which is
why this stage needs the container (ROOT + FCCAnalyses) while the drawing stage
(src/plotting/event_display_plots.py) is plain matplotlib.

Usage (inside the container, see scripts/make_event_displays.sh):
    python src/event_displays.py --max-files 2 --n-per-window 20
"""
import argparse
import os
import pickle
import sys
import types

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

DEFAULT_PROCESS = "p8_ee_ZH_6jet_LF_ecm240"

# Measured on $PATH_TO_HISTOGRAMS/PF_Durham/p8_ee_ZH_6jet_LF_ecm240.root:h_mH_gen
# (mode 124.9 GeV, median 125.9, q05 114.6, q95 162.2). The distribution is
# strongly asymmetric: ~1% of entries below 108 GeV but ~5% above 160 GeV.
DEFAULT_WINDOWS = [
    ("far_left", 90.0, 108.0),
    ("peak", 123.0, 127.0),
    ("far_right", 160.0, 200.0),
]

# Two files of the 6-jet LF sample are anomalous: 0.root has 38637 entries and no
# podio_metadata, 2.root is 636 MB, against ~107 MB / 5000 entries for every
# other file. Skip them so a small --max-files run is representative and fast.
ANOMALOUS_FILES = {"0.root", "2.root"}

# Relative tolerance for the constituent-momentum-sum check (assertion 3).
MOMENTUM_SUM_RTOL = 1e-3


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", default=os.environ.get("PATH_TO_DATASET"),
                   help="Dataset root (default $PATH_TO_DATASET)")
    p.add_argument("--process", default=DEFAULT_PROCESS)
    p.add_argument("--histograms", default=os.environ.get("PATH_TO_HISTOGRAMS"),
                   help="Histogram tree, used only by --windows auto")
    p.add_argument("--method", default="PF_Durham",
                   help="Method directory --windows auto reads the mH histogram from")
    p.add_argument("--hist", default="h_mH_gen",
                   help="Histogram --windows auto derives the windows from")
    p.add_argument("--jet-algorithm", default="Durham",
                   choices=["Durham", "AK", "EEAK", "EECA", "EEKT"])
    p.add_argument("--AK-radius", type=float, default=-1.0)
    p.add_argument("--jet-matching-radius", type=float, default=0.3)
    p.add_argument("--windows", choices=["default", "auto"], default="default")
    p.add_argument("--window", nargs=3, action="append", metavar=("LABEL", "LOW", "HIGH"),
                   help="Override the windows; repeatable")
    p.add_argument("--n-per-window", type=int, default=20)
    p.add_argument("--select", choices=["stratified", "first", "random"],
                   default="stratified",
                   help="Which candidates to keep: stratified spreads them evenly "
                        "over the window instead of piling up at its dense edge")
    p.add_argument("--seed", type=int, default=0, help="only used by --select random")
    p.add_argument("--max-files", type=int, default=2)
    p.add_argument("--start-file", type=int, default=0)
    p.add_argument("--file-list", nargs="+", default=None,
                   help="Explicit ROOT files, overriding --max-files/--start-file")
    p.add_argument("--no-match-histogram-selection", action="store_true",
                   help="Skip the fully-matched-jets filter. Off by default, "
                        "because the PF_Durham tree behind the decomposition "
                        "figures applies it")
    p.add_argument("--output", default=None)
    p.add_argument("--draw", action="store_true",
                   help="Also run the drawing stage in-process")
    return p.parse_args(argv)


def resolve_windows(args):
    """(label, low, high) triples, plus a printed comparison of both recipes."""
    if args.window:
        return [(lab, float(lo), float(hi)) for lab, lo, hi in args.window]
    if args.windows == "default":
        return list(DEFAULT_WINDOWS)

    import uproot

    path = os.path.join(args.histograms, args.method, f"{args.process}.root")
    hist = uproot.open(path)[args.hist]
    y, edges = hist.values(), hist.axis().edges()
    centres = 0.5 * (edges[:-1] + edges[1:])
    total = y.sum()
    if total <= 0:
        raise SystemExit(f"{args.hist} in {path} is empty")
    # Mode on a 2 GeV rebin: the native 0.25 GeV bins are too noisy for an argmax.
    group = max(1, int(round(2.0 / (edges[1] - edges[0]))))
    n_whole = (len(y) // group) * group
    rebinned = y[:n_whole].reshape(-1, group).sum(axis=1)
    rebinned_centres = centres[:n_whole].reshape(-1, group).mean(axis=1)
    mode = float(rebinned_centres[np.argmax(rebinned)])
    cdf = np.cumsum(y) / total

    def q(frac):
        return float(np.interp(frac, cdf, centres))

    windows = [
        ("far_left", q(0.002), q(0.02)),
        ("peak", mode - 2.0, mode + 2.0),
        ("far_right", q(0.95), q(0.995)),
    ]
    print(f"--windows auto, from {path}:{args.hist}")
    print(f"  mode (2 GeV bins) = {mode:.2f} GeV, median = {q(0.5):.2f} GeV")
    for (lab, lo, hi), (_, dlo, dhi) in zip(windows, DEFAULT_WINDOWS):
        print(f"  {lab:10s} auto [{lo:7.2f}, {hi:7.2f}]   hardcoded default "
              f"[{dlo:7.2f}, {dhi:7.2f}]")
    return windows


def select_files(args):
    directory = os.path.join(args.input, args.process)
    if not os.path.isdir(directory):
        raise SystemExit(f"no such process directory: {directory}")
    if args.file_list:
        return list(args.file_list)
    names = sorted(
        (n for n in os.listdir(directory)
         if n.endswith(".root") and n not in ANOMALOUS_FILES),
        key=lambda n: int(n[:-5]) if n[:-5].isdigit() else 0,
    )
    chosen = names[args.start_file: args.start_file + args.max_files]
    if not chosen:
        raise SystemExit(f"no usable ROOT files in {directory}")
    skipped = sorted(ANOMALOUS_FILES & set(os.listdir(directory)))
    if skipped:
        print(f"Skipping known-anomalous input file(s): {', '.join(skipped)}")
    return [os.path.join(directory, n) for n in chosen]


def setup_interpreter():
    """Replicate what `fccanalysis run` does before building the graph."""
    import ROOT

    ROOT.gROOT.SetBatch(True)
    if ROOT.gSystem.Load("libFCCAnalyses") < 0:
        raise SystemExit("could not load libFCCAnalyses - is key4hep sourced?")
    ROOT.dummyLoader()
    ROOT.gInterpreter.Declare("using namespace FCCAnalyses;")
    ROOT.gInterpreter.ProcessLine(".O3")
    for rel in ("histmaker_functions/jet_tools.h",
                "histmaker_functions/utils.h",
                "histmaker_functions/functions.h",
                "histmaker_functions/event_display_functions.h"):
        ROOT.gInterpreter.Declare(f'#include "{os.path.join(REPO, "src", rel)}"')
    # Single-threaded on purpose: Range() and the ordering of the per-window
    # selection are only well defined without implicit MT.
    return ROOT


def to_numpy(value):
    """RVec / std::vector -> numpy array."""
    return np.asarray([v for v in value])


def check_event(event, index, n_jets, label):
    """The consistency checks that make the display trustworthy. Raise, not warn.

    Assertion 3 is the important one: it proves both that a FastJet constituent
    index refers to the input particle collection (JetClustering.cc is not
    shipped on cvmfs, so this cannot be read off the source) and that
    constituent_jet_index reproduces fastjet_to_vec_rp_jet's pT sort. If it ever
    fails, particles would be silently attributed to the wrong jet.
    """
    where = f"{label} event #{index}"
    assign = event["part_jet_index"]
    n_kept = min(n_jets, len(event["jet_eta"]))

    unassigned = int(np.sum(assign < 0))
    assert unassigned == 0, (
        f"{where}: {unassigned} of {len(assign)} particles are in no kept jet; "
        "expected 0 for exclusive Durham clustering to N jets"
    )
    assert set(np.unique(assign)) == set(range(n_kept)), (
        f"{where}: jet indices {sorted(set(np.unique(assign)))} != 0..{n_kept - 1}"
    )

    jet_pt = event["jet_pt"]
    assert np.all(np.diff(jet_pt) <= 1e-4), (
        f"{where}: gen jets are not pT-ordered ({jet_pt}); "
        "fastjet_to_vec_rp_jet's sort assumption is broken"
    )

    for jet in range(n_kept):
        mask = assign == jet
        for axis in ("px", "py", "pz"):
            summed = float(np.sum(event[f"part_{axis}"][mask]))
            expected = float(event[f"jet_{axis}"][jet])
            scale = max(abs(expected), 1.0)
            assert abs(summed - expected) / scale < MOMENTUM_SUM_RTOL, (
                f"{where}: jet {jet} {axis} from constituents = {summed:.4f} but "
                f"the jet has {expected:.4f}. The constituent->particle index "
                "convention or the pT sort does not hold; the display would "
                "attribute particles to the wrong jets."
            )

    higgs = [j for j in event["parton_to_jet"] if j >= 0]
    if len(higgs) == len(event["parton_to_jet"]):
        e = sum(float(event["jet_energy"][j]) for j in higgs)
        px, py, pz = (sum(float(event[f"jet_{a}"][j]) for j in higgs)
                      for a in ("px", "py", "pz"))
        mass2 = e * e - (px * px + py * py + pz * pz)
        mass = np.sqrt(mass2) if mass2 > 0 else 0.0
        assert abs(mass - event["mH_gen"]) < 1e-2, (
            f"{where}: mH recomputed from the Higgs-matched jets = {mass:.4f} but "
            f"inv_mass_gen = {event['mH_gen']:.4f}"
        )


def pick(indices, values, n, how, seed):
    """Which candidates to keep, deterministically unless how == 'random'."""
    if len(indices) <= n:
        return list(range(len(indices)))
    if how == "first":
        return list(range(n))
    if how == "random":
        rng = np.random.default_rng(seed)
        return sorted(rng.choice(len(indices), size=n, replace=False).tolist())
    order = np.argsort(values)
    return sorted(order[np.linspace(0, len(order) - 1, n).round().astype(int)].tolist())


def main(argv=None):
    args = parse_args(argv)
    if not args.input:
        raise SystemExit("set PATH_TO_DATASET or pass --input")

    from src.process_config import NUMBER_OF_JETS, NUMBER_OF_HIGGS_JETS
    n_jets = NUMBER_OF_JETS[args.process]
    n_higgs_jets = NUMBER_OF_HIGGS_JETS[args.process]

    windows = resolve_windows(args)
    files = select_files(args)
    print(f"Process {args.process}: {n_jets} jets, {n_higgs_jets} Higgs partons")
    print(f"Reading {len(files)} file(s): {', '.join(os.path.basename(f) for f in files)}")

    ROOT = setup_interpreter()
    from src.histmaker_tools.event_display import (
        build_event_display_graph, PAYLOAD_COLUMNS, SCALAR_COLUMNS,
    )

    graph_args = types.SimpleNamespace(
        jet_algorithm=args.jet_algorithm,
        AK_radius=args.AK_radius,
        jet_matching_radius=args.jet_matching_radius,
        energy_recovery=False,
        ideal_matching=False,
    )
    chain = ROOT.std.vector("string")()
    for path in files:
        chain.push_back(path)
    df = ROOT.RDataFrame("events", chain)
    n_input = df.Count()
    df = build_event_display_graph(
        df, graph_args, n_jets, n_higgs_jets,
        apply_matched_filter=not args.no_match_histogram_selection,
    )
    # mH is -1 when the Higgs partons did not all land on distinct jets
    # (functions.h invariant_mass), which is exactly the "undefined" case.
    df = df.Filter("mH_gen > 0", "mH defined")

    selections = {}
    for label, low, high in windows:
        node = df.Filter(f"mH_gen >= {low} && mH_gen < {high}", label)
        selections[label] = (node.Count(), node)

    total_input = n_input.GetValue()
    print(f"\n{total_input} events read")
    report = df.Report()
    report.Print()

    payload = {"process": args.process, "windows": [], "n_input_events": total_input,
               "jet_algorithm": args.jet_algorithm, "n_jets": n_jets,
               "n_higgs_jets": n_higgs_jets, "hist": args.hist,
               "matched_filter": not args.no_match_histogram_selection}

    for label, low, high in windows:
        count_ptr, node = selections[label]
        n_candidates = count_ptr.GetValue()
        print(f"\n--- window {label} [{low}, {high}] GeV: {n_candidates} candidate(s)")
        if n_candidates == 0:
            print("    no events in this window - increase --max-files")
            payload["windows"].append(
                {"label": label, "low": low, "high": high, "events": [],
                 "n_candidates": 0})
            continue
        columns = node.AsNumpy(PAYLOAD_COLUMNS)
        events = []
        for i in range(n_candidates):
            event = {}
            for name in PAYLOAD_COLUMNS:
                value = columns[name][i]
                event[name] = value if name in SCALAR_COLUMNS else to_numpy(value)
            events.append(event)
        keep = pick(list(range(len(events))),
                    [float(e["mH_gen"]) for e in events],
                    args.n_per_window, args.select, args.seed)
        kept = [events[i] for i in keep]
        for i, event in enumerate(kept):
            assert low <= event["mH_gen"] < high, (
                f"{label} event #{i}: mH {event['mH_gen']} outside [{low}, {high})"
            )
            check_event(event, i, n_jets, label)
        rate = n_candidates / total_input * 100
        print(f"    keeping {len(kept)} ({args.select}); "
              f"yield {rate:.2f}% of the {total_input} events read")
        payload["windows"].append(
            {"label": label, "low": low, "high": high, "events": kept,
             "n_candidates": n_candidates, "yield_percent": rate})

    output = args.output or os.path.join(
        os.environ.get("PATH_TO_HISTOGRAMS", "."), "plots", "event_displays",
        f"{args.process}_payload.pkl")
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "wb") as handle:
        pickle.dump(payload, handle)
    n_total = sum(len(w["events"]) for w in payload["windows"])
    print(f"\nAll consistency checks passed. Wrote {n_total} event(s) to {output}")

    if args.draw:
        from src.plotting.event_display_plots import draw_payload
        pdf = os.path.splitext(output)[0] + ".pdf"
        draw_payload(payload, pdf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
