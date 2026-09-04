"""Trial 15 - three_param weighted with a 1% relative error FLOOR.

sigma_i = sqrt(stat_i^2 + (0.01*y_i)^2). The floor stops the millions-of-events
mid-E bins from utterly dominating the weight, so low/high-E points still
constrain the fit. Should bring chi2/ndf down toward O(1) and balance the fit.

RESULTS: MAE=8.9311e-04 RMSE=1.3878e-03 maxAE=3.6781e-03 looCV_MAE=9.9555e-04 chi2/ndf=3.058e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 3:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.01)
    try:
        popt, _ = curve_fit(m_stoch_const_noise, x, y, p0=[0.5, 0.03, 0.1],
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("15_floor_1pct",
                      "three_param weighted, 1% relative error floor", fit_one)
