"""Trial 09 - four-param log(E) model, error-weighted.

Trial 06's model but weighted by stored stat errors. Checks whether the extra
log DOF plus weighting improves the (statistically meaningful) chi2 without
wrecking point-wise error.

RESULTS: MAE=6.7508e-04 RMSE=1.2410e-03 maxAE=3.4065e-03 looCV_MAE=6.5458e-04 chi2/ndf=7.121e+03 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_logE


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 4:
        return None
    sigma = harness.safe_sigma(err, y, floor_frac=0.0)
    try:
        popt, _ = curve_fit(m_logE, x, y, p0=[0.5, 0.03, 0.1, 0.0],
                            sigma=sigma, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    return lambda X: m_logE(X, *popt)


if __name__ == "__main__":
    harness.run_trial("09_logE_4param_weighted",
                      "a/sqrt(E)+b+c/E+d*log(E), stat-error weighted", fit_one)
