"""Trial 35 - trial-20 model (e/E^2) via robust soft_l1 loss.

The winning 4-param model, but fit with an outlier-resistant soft_l1 loss so a
single noisy turn-up/tail bin can't dominate. Tests whether robustness shaves
the max abs error (trial 20's weak spot vs its great MAE).

RESULTS: MAE=3.9060e-04 RMSE=4.9529e-04 maxAE=1.0197e-03 looCV_MAE=8.2026e-04 chi2/ndf=1.174e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None

    def resid(p):
        return m_stoch_const_noise2(x, *p) - y

    try:
        res = least_squares(resid, x0=[0.5, 0.03, 0.1, 0.01], loss="soft_l1",
                            f_scale=np.median(np.abs(y)) or 1.0, max_nfev=40000)
    except Exception:
        return None
    if not res.success and res.status <= 0:
        return None
    return lambda X: m_stoch_const_noise2(X, *res.x)


if __name__ == "__main__":
    harness.run_trial("35_invE2_soft_l1",
                      "a/sqrt(E)+b+c/E+e/E^2, robust soft_l1", fit_one)
