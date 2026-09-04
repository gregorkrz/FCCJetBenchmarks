"""Trial 18 - three_param via robust soft_l1 loss (outlier-resistant).

Uses least_squares with loss='soft_l1' so a single bad point (e.g. a low-stat
turn-up bin) can't dominate. Unweighted residuals.

RESULTS: MAE=8.3797e-04 RMSE=1.0399e-03 maxAE=2.2245e-03 looCV_MAE=1.1694e-03 chi2/ndf=5.488e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None

    def resid(p):
        return m_stoch_const_noise(x, *p) - y

    try:
        res = least_squares(resid, x0=[0.5, 0.03, 0.1], loss="soft_l1",
                            f_scale=np.median(np.abs(y)) or 1.0, max_nfev=20000)
    except Exception:
        return None
    if not res.success and res.status <= 0:
        return None
    popt = res.x
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("18_robust_soft_l1",
                      "three_param, robust soft_l1 loss", fit_one)
