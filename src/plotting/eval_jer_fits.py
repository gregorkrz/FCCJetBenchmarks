"""Quick offline evaluation of the energy-JER fit methods from dashboard_data.json.

Reads the packaged dashboard JSON (no ROOT, no refitting), evaluates every
stored fit method at the actual sigma/E data points, and reports goodness-of-fit
metrics (MAE, MSE/RMSE, max abs error, and - where per-point errors are stored -
a reduced chi2) so the fitting methods can be compared head to head.

Usage:
    python src/plotting/eval_jer_fits.py \
        --json $PATH_TO_HISTOGRAMS/plots/dashboard_data.json
    # optionally restrict:  --method PF_Durham --process p8_ee_ZH_bbbb_ecm240
"""
import argparse
import json

import numpy as np


# popt is in curve_fit PARAMETER order, i.e. exactly the order of the model
# function's args (NOT the A/B/C display labels, which swap the last two):
#   three_param = _model_three_param(E, a, b, c) = a/sqrt(E) + b + c/E
#       -> popt = [a (stochastic), b (constant), c (noise)]
#   two_param   = _model_two_param(E, a, b)      = a/sqrt(E) + b
#       -> popt = [a (stochastic), b (constant)]
# The "_weighted" variants share these forms.
def _model_value(base_model, E, popt):
    E = np.asarray(E, dtype=float)
    if base_model == "three_param":
        a, b, c = popt
        return a / np.sqrt(E) + b + c / E
    if base_model == "two_param":
        a, b = popt
        return a / np.sqrt(E) + b
    raise ValueError(f"Unknown base model: {base_model}")


def _base_of(fit_name):
    return fit_name.replace("_weighted", "")


def evaluate_entry(entry):
    """Return {fit_name: metrics} for one energy._all entry."""
    x = np.asarray(entry.get("mid_points", []), dtype=float)
    y = np.asarray(entry.get("sigma_over_E", []), dtype=float)
    err = entry.get("sigma_over_E_err")
    err = np.asarray(err, dtype=float) if err is not None else None
    out = {}
    if x.size < 2:
        return out
    for fit_name, fit in (entry.get("fits") or {}).items():
        if not fit or "popt" not in fit:
            continue
        try:
            pred = _model_value(_base_of(fit_name), x, fit["popt"])
        except Exception:
            continue
        resid = y - pred
        n = resid.size
        metrics = {
            "MAE": float(np.mean(np.abs(resid))),
            "MSE": float(np.mean(resid ** 2)),
            "RMSE": float(np.sqrt(np.mean(resid ** 2))),
            "maxAE": float(np.max(np.abs(resid))),
            "n": int(n),
        }
        if err is not None and err.shape == resid.shape:
            good = np.isfinite(err) & (err > 0)
            ndf = int(good.sum()) - len(fit["popt"])
            if ndf > 0:
                metrics["chi2_ndf"] = float(np.sum((resid[good] / err[good]) ** 2) / ndf)
        out[fit_name] = metrics
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True, help="Path to dashboard_data.json")
    ap.add_argument("--method", default=None, help="Restrict to one method")
    ap.add_argument("--process", default=None, help="Restrict to one process")
    ap.add_argument("--metric", default="MAE",
                    choices=["MAE", "MSE", "RMSE", "maxAE", "chi2_ndf"],
                    help="Metric used for the per-method aggregate ranking")
    args = ap.parse_args()

    with open(args.json) as fh:
        d = json.load(fh)

    # Aggregate metric values per fit method across all (method, process) entries.
    agg = {}
    rows = []
    for method, mv in d["methods"].items():
        if args.method and method != args.method:
            continue
        for process, pv in mv["processes"].items():
            if args.process and process != args.process:
                continue
            entry = (pv.get("energy") or {}).get("_all")
            if not entry:
                continue
            res = evaluate_entry(entry)
            for fit_name, metrics in res.items():
                agg.setdefault(fit_name, []).append(metrics)
            if args.method or args.process:
                rows.append((method, process, res))

    # Per-(method,process) detail, only when a filter is applied (else too noisy).
    for method, process, res in rows:
        print(f"\n{method} / {process}")
        for fit_name, m in sorted(res.items()):
            extra = f"  chi2/ndf={m['chi2_ndf']:.3g}" if "chi2_ndf" in m else ""
            print(f"  {fit_name:24s}  MAE={m['MAE']:.3e}  RMSE={m['RMSE']:.3e}  "
                  f"maxAE={m['maxAE']:.3e}{extra}")

    # Aggregate summary across everything selected.
    print("\n=== aggregate across "
          f"{sum(len(v) for v in agg.values()) // max(len(agg), 1)} entries "
          f"(mean over entries), ranked by mean {args.metric} ===")
    summary = []
    for fit_name, ms in agg.items():
        row = {"fit": fit_name, "count": len(ms)}
        for key in ["MAE", "MSE", "RMSE", "maxAE"]:
            row[key] = float(np.mean([m[key] for m in ms]))
        chi = [m["chi2_ndf"] for m in ms if "chi2_ndf" in m]
        row["chi2_ndf"] = float(np.mean(chi)) if chi else None
        summary.append(row)

    def sort_key(r):
        v = r.get(args.metric)
        return (v is None, v if v is not None else 0.0)

    summary.sort(key=sort_key)
    hdr = f"{'method':24s} {'count':>6s} {'MAE':>12s} {'RMSE':>12s} {'maxAE':>12s} {'chi2/ndf':>12s}"
    print(hdr)
    print("-" * len(hdr))
    for r in summary:
        chi = f"{r['chi2_ndf']:.4g}" if r["chi2_ndf"] is not None else "-"
        print(f"{r['fit']:24s} {r['count']:6d} {r['MAE']:12.4e} {r['RMSE']:12.4e} "
              f"{r['maxAE']:12.4e} {chi:>12s}")


if __name__ == "__main__":
    main()
