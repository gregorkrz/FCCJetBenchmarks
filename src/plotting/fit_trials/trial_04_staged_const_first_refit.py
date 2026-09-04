"""Trial 04 - STAGED then RELEASE: constant-first seed, then full refit.

Like trial 03 but the fixed-b stage is only used to seed a good starting point;
step 3 releases all three params and refits a/sqrt(E)+b+c/E from that seed. Tests
whether the staged approach helps mainly via better initialisation.

RESULTS: MAE=8.4170e-04 RMSE=1.0386e-03 maxAE=2.1871e-03 looCV_MAE=1.1685e-03 chi2/ndf=5.753e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise

N_TAIL = 3


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    b0 = float(np.mean(ys[-N_TAIL:]))
    try:
        p_ac, _ = curve_fit(lambda E, a, c: a / np.sqrt(E) + b0 + c / E,
                            xs, ys, p0=[0.5, 0.1], maxfev=20000)
        a0, c0 = p_ac
        popt, _ = curve_fit(m_stoch_const_noise, xs, ys, p0=[a0, b0, c0], maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("04_staged_const_first_refit",
                      f"seed b=mean(last {N_TAIL}), fit a,c, then release+refit all 3", fit_one)
