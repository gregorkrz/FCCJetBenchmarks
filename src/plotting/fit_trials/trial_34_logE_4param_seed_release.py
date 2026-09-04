"""Trial 34 - trial-06 log model with staged-b seed then release. 4 params.

Trial 06 (a/sqrt(E)+b+c/E+d*log(E)) was round 1's #2. Give it the same staged-b
seeding that helped elsewhere: pin b from the tail to seed a,c,d, then release
all four. Checks whether the log model was starting-point limited.

RESULTS: MAE=4.5401e-04 RMSE=5.6089e-04 maxAE=1.1059e-03 looCV_MAE=8.3772e-04 chi2/ndf=1.609e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE

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
            lambda E, a, c, d: a / np.sqrt(E) + b0 + c / E + d * np.log(E),
            xs, ys, p0=[0.5, 0.1, 0.0], maxfev=40000)
        a0, c0, d0 = p_seed
        popt, _ = curve_fit(m_logE, xs, ys, p0=[a0, b0, c0, d0], maxfev=40000)
    except Exception:
        return None
    return lambda X: m_logE(X, *popt)


if __name__ == "__main__":
    harness.run_trial("34_logE_4param_seed_release",
                      "log 4-param, staged-b seed then release", fit_one)
