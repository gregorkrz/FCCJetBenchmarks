"""Trial 41 - the 5-param log+e/E^2 champion, error-weighted with 1% floor.

Trial 26 (current MAE leader) was unweighted. Add balanced 1% weighting to see
if the statistically-meaningful chi2 drops without hurting the point-wise error.

RESULTS: MAE=2.7330e-04 RMSE=3.5670e-04 maxAE=7.9870e-04 looCV_MAE=9.6101e-04 chi2/ndf=6.074e+03 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE_and_invE2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 6:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.01)
    try:
        popt, _ = curve_fit(m_logE_and_invE2, x, y, p0=[0.5, 0.03, 0.1, 0.01, 0.0],
                            sigma=sigma, absolute_sigma=True, maxfev=60000)
    except Exception:
        return None
    return lambda X: m_logE_and_invE2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("41_logE_invE2_weighted",
                      "5-param log+e/E^2, weighted 1% floor", fit_one)
