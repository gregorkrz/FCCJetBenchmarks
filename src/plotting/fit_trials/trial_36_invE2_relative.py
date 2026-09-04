"""Trial 36 - trial-20 model (e/E^2), minimising RELATIVE residuals (sigma=y).

The winning model fit for uniform FRACTIONAL accuracy across the energy range
(weight 1/y). Useful if what matters is % error on sigma/E at every E, not
absolute error (which favours the big low-E points).

RESULTS: MAE=4.4871e-04 RMSE=6.0093e-04 maxAE=1.3363e-03 looCV_MAE=8.7901e-04 chi2/ndf=1.056e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    sigma = np.abs(y) if np.all(y > 0) else None
    try:
        popt, _ = curve_fit(m_stoch_const_noise2, x, y, p0=[0.5, 0.03, 0.1, 0.01],
                            sigma=sigma, absolute_sigma=False, maxfev=40000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("36_invE2_relative",
                      "a/sqrt(E)+b+c/E+e/E^2, relative residuals", fit_one)
