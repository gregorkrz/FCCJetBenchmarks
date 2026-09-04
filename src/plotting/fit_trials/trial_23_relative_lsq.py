"""Trial 23 - three_param minimising RELATIVE residuals (resid/y).

Weights every point by 1/y, i.e. minimises fractional error. Natural when you
care about % accuracy of sigma/E uniformly across the energy range rather than
absolute error (which favours the large low-E values).

RESULTS: MAE=8.9263e-04 RMSE=1.3825e-03 maxAE=3.6608e-03 looCV_MAE=9.9837e-04 chi2/ndf=3.076e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3 or np.any(y <= 0):
        # fall back to plain LSQ if any non-positive y
        sigma = None
    else:
        sigma = np.abs(y)
    try:
        popt, _ = curve_fit(m_stoch_const_noise, x, y, p0=[0.5, 0.03, 0.1],
                            sigma=sigma, absolute_sigma=False, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("23_relative_lsq",
                      "three_param, relative residuals (sigma=y)", fit_one)
