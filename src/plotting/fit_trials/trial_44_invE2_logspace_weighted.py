"""Trial 44 - trial-20 model (e/E^2) in log space, stat-error weighted.

Log-space fit (trial 38) but weighting each log-residual by the relative error
delta(y)/y ~ 1/sqrt(2N), i.e. proper error propagation into log space. The most
statistically principled version of the strong 4-param e/E^2 model.

RESULTS: MAE=4.4772e-04 RMSE=5.8928e-04 maxAE=1.2839e-03 looCV_MAE=8.8258e-04 chi2/ndf=1.051e+04 (fit 159, failed 0)
"""
import os, sys
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness
from models import m_stoch_const_noise2


def fit_one(x, y, err):
    x, y, err = harness.clean(x, y, err)
    if x.size < 5 or np.any(y <= 0):
        return None
    # delta(log y) = delta(y)/y; add a small floor so high-N points don't dominate.
    if err is not None:
        rel = np.where(np.isfinite(err) & (err > 0), err / np.abs(y), np.nan)
        sigma_log = np.sqrt(np.nan_to_num(rel, nan=0.0) ** 2 + 0.01 ** 2)
    else:
        sigma_log = None

    def log_model(E, a, b, c, e):
        return np.log(m_stoch_const_noise2(E, a, b, c, e))

    try:
        popt, _ = curve_fit(log_model, x, np.log(y), p0=[0.5, 0.03, 0.1, 0.01],
                            sigma=sigma_log, absolute_sigma=True, maxfev=40000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise2(X, *popt)


if __name__ == "__main__":
    harness.run_trial("44_invE2_logspace_weighted",
                      "e/E^2 model, log space, error-weighted", fit_one)
