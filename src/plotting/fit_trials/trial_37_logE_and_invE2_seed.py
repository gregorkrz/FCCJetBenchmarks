"""Trial 37 - the 5-param log+e/E^2 model with staged-b seed then release.

Trial 26's 5-param kitchen-sink model, but seeded via a fixed-b sub-fit first to
help the extra params converge. The most expressive model tried, given the best
initialisation strategy found so far.

RESULTS: MAE=2.5478e-04 RMSE=3.2658e-04 maxAE=6.7990e-04 looCV_MAE=8.9428e-04 chi2/ndf=6.454e+03 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE_and_invE2

N_TAIL = 3


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 6:
        return None
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    b0 = float(np.mean(ys[-N_TAIL:]))
    try:
        p_seed, _ = curve_fit(
            lambda E, a, c, e, d: a / np.sqrt(E) + b0 + c / E + e / E ** 2 + d * np.log(E),
            xs, ys, p0=[0.5, 0.1, 0.01, 0.0], maxfev=60000)
        a0, c0, e0, d0 = p_seed
        popt, _ = curve_fit(m_logE_and_invE2, xs, ys, p0=[a0, b0, c0, e0, d0], maxfev=60000)
    except Exception:
        return None
    return lambda X: m_logE_and_invE2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("37_logE_and_invE2_seed",
                      "5-param log+e/E^2, staged-b seed then release", fit_one)
