"""Trial 42 - 5-param log+e/E^2 champion, relative residuals (sigma=y).

Trial 26's model fit for uniform fractional accuracy. Checks whether the leader
can be pushed further on % error while keeping absolute error competitive.

RESULTS: MAE=2.7325e-04 RMSE=3.5633e-04 maxAE=7.9678e-04 looCV_MAE=9.6248e-04 chi2/ndf=6.104e+03 (fit 159, failed 0)
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
    sigma = np.abs(y) if np.all(y > 0) else None
    try:
        popt, _ = curve_fit(m_logE_and_invE2, x, y, p0=[0.5, 0.03, 0.1, 0.01, 0.0],
                            sigma=sigma, absolute_sigma=False, maxfev=60000)
    except Exception:
        return None
    return lambda X: m_logE_and_invE2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("42_logE_and_invE2_relative",
                      "5-param log+e/E^2, relative residuals", fit_one)
