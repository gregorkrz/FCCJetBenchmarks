"""Trial 02 - three_param linear sum, error-weighted (Option A), unbounded.

Same model as trial 01 but weighted by the stored per-point stat error
delta(sigma/E)=(sigma/E)/sqrt(2N). Tests whether pure statistical weighting
(dominated by high-N mid-E bins) helps or hurts point-wise error.

RESULTS: MAE=1.0871e-03 RMSE=2.0752e-03 maxAE=5.9348e-03 looCV_MAE=8.7278e-04 chi2/ndf=1.766e+04 (fit 159, failed 0)
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
    sigma = harness.safe_sigma(err, y, floor_frac=0.0)
    try:
        popt, _ = curve_fit(m_stoch_const_noise, x, y, p0=[0.5, 0.03, 0.1],
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("02_three_param_weighted",
                      "a/sqrt(E)+b+c/E, stat-error weighted (Option A), unbounded", fit_one)
