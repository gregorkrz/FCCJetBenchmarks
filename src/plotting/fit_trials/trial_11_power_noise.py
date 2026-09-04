"""Trial 11 - a/sqrt(E)+b+c*E^(-p)  (free exponent on the noise-like term).

Keeps the stochastic 1/sqrt(E) fixed but lets the steep low-E term have a free
exponent (nominally 1 for noise). 4 params.

RESULTS: MAE=1.3415e-03 RMSE=1.6988e-03 maxAE=3.7364e-03 looCV_MAE=1.5738e-03 chi2/ndf=7.537e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_power_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    try:
        popt, _ = curve_fit(m_power_noise, x, y, p0=[0.5, 0.03, 0.1, 1.0],
                            bounds=([0, 0, 0, 0.3], [np.inf, np.inf, np.inf, 3.0]), maxfev=20000)
    except Exception:
        return None
    return lambda X: m_power_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("11_power_noise",
                      "a/sqrt(E)+b+c*E^(-p), free noise exponent p in [0.3,3]", fit_one)
