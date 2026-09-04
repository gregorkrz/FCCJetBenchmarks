"""Trial 19 - three_param via robust Huber loss.

Like trial 18 but Huber loss (quadratic near 0, linear in the tails) - a milder
robustification than soft_l1.

RESULTS: MAE=8.4095e-04 RMSE=1.0387e-03 maxAE=2.1957e-03 looCV_MAE=1.1715e-03 chi2/ndf=5.679e+04 (fit 159, failed 0)
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
        res = least_squares(resid, x0=[0.5, 0.03, 0.1], loss="huber",
                            f_scale=np.median(np.abs(y)) or 1.0, max_nfev=20000)
    except Exception:
        return None
    if not res.success and res.status <= 0:
        return None
    return lambda X: m_stoch_const_noise(X, *res.x)


if __name__ == "__main__":
    harness.run_trial("19_robust_huber",
                      "three_param, robust Huber loss", fit_one)
