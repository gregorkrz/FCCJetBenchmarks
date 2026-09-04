"""Trial 03 - STAGED: fix constant b from high-E points, then fit a,c.

User's suggestion. Step 1: estimate the constant term b as the mean of sigma/E
over the last few (highest-E) points, where a/sqrt(E) and c/E have died away so
sigma/E ~ b. Step 2: fix b and fit the remaining a/sqrt(E)+c/E to all points.

RESULTS: MAE=2.6909e-03 RMSE=3.2159e-03 maxAE=6.8081e-03 looCV_MAE=3.0265e-03 chi2/ndf=2.827e+05 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness

N_TAIL = 3  # number of highest-E points used to pin the constant


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    order = np.argsort(x)
    x, y = x[order], y[order]
    b = float(np.mean(y[-N_TAIL:]))

    def resid_model(E, a, c):
        return a / np.sqrt(E) + b + c / E

    try:
        popt, _ = curve_fit(resid_model, x, y, p0=[0.5, 0.1], maxfev=20000)
    except Exception:
        return None
    a, c = popt
    return lambda X: a / np.sqrt(np.asarray(X, float)) + b + c / np.asarray(X, float)


if __name__ == "__main__":
    harness.run_trial("03_staged_const_first",
                      f"fix b=mean(last {N_TAIL} pts), then fit a/sqrt(E)+c/E", fit_one)
