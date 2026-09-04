"""Trial 40 - a/sqrt(E)+b+c/E^r (free steep exponent), weighted 1% floor.

Trial 30's free-exponent steep term with balanced weighting. Lets the data pick
r AND respects statistical weight - the most flexible 4-param low-E treatment.

RESULTS: MAE=1.4174e-03 RMSE=2.1640e-03 maxAE=5.2398e-03 looCV_MAE=1.3982e-03 chi2/ndf=5.599e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_free_low_power


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.01)
    try:
        popt, _ = curve_fit(m_free_low_power, x, y, p0=[0.5, 0.03, 0.1, 1.5],
                            bounds=([0, 0, 0, 1.0], [np.inf, np.inf, np.inf, 3.0]),
                            sigma=sigma, absolute_sigma=True, maxfev=40000)
    except Exception:
        return None
    return lambda X: m_free_low_power(X, *popt)


if __name__ == "__main__":
    harness.run_trial("40_free_low_power_weighted",
                      "a/sqrt(E)+b+c/E^r, free r, weighted 1% floor", fit_one)
