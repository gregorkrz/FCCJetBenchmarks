"""Trial 22 - STAGED with b fixed from only the last 2 points, kept fixed.

Variant of trial 03 exploring tail length: use just the 2 highest-E points for
the plateau estimate and never release b. Tests sensitivity of the staged idea
to how many tail points define the constant.

RESULTS: MAE=2.6460e-03 RMSE=3.1648e-03 maxAE=6.7142e-03 looCV_MAE=2.9683e-03 chi2/ndf=2.738e+05 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness

N_TAIL = 2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    order = np.argsort(x)
    x, y = x[order], y[order]
    b = float(np.mean(y[-N_TAIL:]))
    try:
        popt, _ = curve_fit(lambda E, a, c: a / np.sqrt(E) + b + c / E,
                            x, y, p0=[0.5, 0.1], maxfev=20000)
    except Exception:
        return None
    a, c = popt
    return lambda X: a / np.sqrt(np.asarray(X, float)) + b + c / np.asarray(X, float)


if __name__ == "__main__":
    harness.run_trial("22_staged_const_fixed_2pt",
                      f"fix b=mean(last {N_TAIL}), fit a,c (no release)", fit_one)
