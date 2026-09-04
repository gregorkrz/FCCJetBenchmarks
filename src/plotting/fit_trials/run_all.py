"""Run every trial_*.py in this folder and print a ranked leaderboard.

Each trial appends its aggregate metrics to results.jsonl (de-duped by name).
After running them all, this prints a table sorted by the chosen metric and
back-fills the "RESULTS:" line in each trial script's docstring so every script
carries its own headline numbers.

Usage:
    python src/plotting/fit_trials/run_all.py            # rank by RMSE
    python src/plotting/fit_trials/run_all.py --metric MAE
"""
import argparse
import glob
import importlib.util
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(HERE, "results.jsonl")


def _load_module(path):
    spec = importlib.util.spec_from_file_location(
        "trial_" + os.path.basename(path)[:-3], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _backfill_results_line(path, agg):
    with open(path) as fh:
        text = fh.read()
    loo = f"{agg['loo_mae']:.4e}" if agg.get("loo_mae") is not None else "n/a"
    headline = (f"RESULTS: MAE={agg['MAE']:.4e} RMSE={agg['RMSE']:.4e} "
                f"maxAE={agg['maxAE']:.4e} looCV_MAE={loo} "
                f"chi2/ndf={agg['chi2_ndf']:.3e} "
                f"(fit {agg['n_series']}, failed {agg['n_failed']})")
    new = re.sub(r"RESULTS:.*", headline, text, count=1)
    if new != text:
        with open(path, "w") as fh:
            fh.write(new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", default="RMSE",
                    choices=["MAE", "MSE", "RMSE", "maxAE", "chi2_ndf", "loo_mae"])
    args = ap.parse_args()

    trial_paths = sorted(glob.glob(os.path.join(HERE, "trial_*.py")))
    print(f"Running {len(trial_paths)} trials...\n")
    results = []
    for path in trial_paths:
        mod = _load_module(path)
        if not hasattr(mod, "fit_one"):
            continue
        # Re-run via run_trial to persist; grab the returned agg.
        name = None
        desc = ""
        # Trials call run_trial in __main__; call it here directly instead.
        # Recover name/desc from the module by a light convention: reuse the
        # file's own __main__ block by executing run_trial with introspected args
        # is overkill - instead just call harness.run_trial through the module.
        import harness
        # Parse the run_trial(...) call args from the source for name/desc.
        with open(path) as fh:
            src = fh.read()
        m = re.search(r'run_trial\(\s*["\'](.+?)["\']\s*,\s*["\'](.+?)["\']', src, re.S)
        if m:
            name, desc = m.group(1), m.group(2)
        else:
            name = os.path.basename(path)[:-3]
        agg = harness.run_trial(name, desc, mod.fit_one, verbose=False)
        _backfill_results_line(path, agg)
        results.append(agg)
        print(f"  ok  {name:34s} RMSE={agg['RMSE']:.4e} MAE={agg['MAE']:.4e}")

    # Leaderboard
    def key(r):
        v = r.get(args.metric)
        return (v is None, v if v is not None else 0.0)

    results.sort(key=key)
    print(f"\n=== LEADERBOARD (ranked by {args.metric}, lower is better) ===")
    hdr = (f"{'#':>2} {'trial':34s} {'MAE':>11s} {'RMSE':>11s} {'maxAE':>11s} "
           f"{'looCV_MAE':>11s} {'chi2/ndf':>11s} {'fail':>4s}")
    print(hdr)
    print("-" * len(hdr))
    for i, r in enumerate(results, 1):
        chi = f"{r['chi2_ndf']:.3e}" if r.get("chi2_ndf") is not None else "-"
        loo = f"{r['loo_mae']:.4e}" if r.get("loo_mae") is not None else "-"
        print(f"{i:2d} {r['trial']:34s} {r['MAE']:11.4e} {r['RMSE']:11.4e} "
              f"{r['maxAE']:11.4e} {loo:>11s} {chi:>11s} {r['n_failed']:4d}")


if __name__ == "__main__":
    main()
