"""Trial 43 - the no-sqrt 4-param model (trial 33) with staged-b seed.

Trial 33 (b+c/E+e/E^2+d*log(E), no stochastic term) did nearly as well as the
5-param leader with one fewer param - a parsimony win if CV agrees. Give it the
staged-b seed for robust convergence and see if it holds up under LOO CV.

RESULTS: MAE=4.0570e-04 RMSE=5.1114e-04 maxAE=1.0371e-03 looCV_MAE=8.2960e-04 chi2/ndf=1.243e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE_invE2_noStoch

N_TAIL = 3


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    b0 = float(np.mean(ys[-N_TAIL:]))
    try:
        p_seed, _ = curve_fit(
            lambda E, c, e, d: b0 + c / E + e / E ** 2 + d * np.log(E),
            xs, ys, p0=[0.5, 0.05, 0.0], maxfev=40000)
        c0, e0, d0 = p_seed
        popt, _ = curve_fit(m_logE_invE2_noStoch, xs, ys, p0=[b0, c0, e0, d0], maxfev=40000)
    except Exception:
        return None
    return lambda X: m_logE_invE2_noStoch(X, *popt)


if __name__ == "__main__":
    harness.run_trial("43_logE_invE2_noStoch_seed",
                      "b+c/E+e/E^2+d*log(E), staged-b seed (4-param, no sqrt)", fit_one)
