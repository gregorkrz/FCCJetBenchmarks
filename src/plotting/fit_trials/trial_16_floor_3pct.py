"""Trial 16 - three_param weighted with a 3% relative error FLOOR.

Same idea as trial 15 with a larger floor - closer to unweighted behaviour, a
scan point to see how the floor size trades chi2 against point-wise error.

RESULTS: MAE=8.9269e-04 RMSE=1.3831e-03 maxAE=3.6628e-03 looCV_MAE=9.9805e-04 chi2/ndf=3.074e+04 (fit 159, failed 0)
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
    sigma = harness.safe_sigma(err, y, floor_frac=0.03)
    try:
        popt, _ = curve_fit(m_stoch_const_noise, x, y, p0=[0.5, 0.03, 0.1],
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("16_floor_3pct",
                      "three_param weighted, 3% relative error floor", fit_one)
