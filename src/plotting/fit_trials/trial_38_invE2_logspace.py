"""Trial 38 - trial-20 model fit in LOG space: minimise resid of log(sigma/E).

Fitting log(y) turns multiplicative (fractional) error into additive error and
naturally balances the 20x dynamic range of sigma/E across energy - a common
trick for resolution curves. Model predictions are exp'd back for metric eval.

RESULTS: MAE=4.4726e-04 RMSE=5.8685e-04 maxAE=1.2747e-03 looCV_MAE=8.8443e-04 chi2/ndf=1.057e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5 or np.any(y <= 0):
        return None

    def log_model(E, a, b, c, e):
        return np.log(m_stoch_const_noise2(E, a, b, c, e))

    try:
        popt, _ = curve_fit(log_model, x, np.log(y), p0=[0.5, 0.03, 0.1, 0.01], maxfev=40000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("38_invE2_logspace",
                      "a/sqrt(E)+b+c/E+e/E^2, fit in log space", fit_one)
