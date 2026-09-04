"""Trial 45 - the CURRENT PRODUCTION fit, reproduced for a fair reference row.

Exactly what resolution_plots.py / build_dashboard_data.py do for the unweighted
three_param: a/sqrt(E)+b+c/E with bounds b in [0.005, 0.04] and a,c >= 0. This is
the incumbent every other trial should be measured against, on the same footing
and with the same LOO CV as the experiments.

RESULTS: MAE=1.8395e-03 RMSE=2.2906e-03 maxAE=5.1642e-03 looCV_MAE=2.0286e-03 chi2/ndf=1.636e+05 (fit 159, failed 0)
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
    try:
        popt, _ = curve_fit(m_stoch_const_noise, x, y, p0=[0.5, 0.03, 0.1],
                            bounds=([0.0, 0.005, 0.0], [np.inf, 0.04, np.inf]), maxfev=20000)
    except Exception:
        return None
    return lambda X: m_stoch_const_noise(X, *popt)


if __name__ == "__main__":
    harness.run_trial("45_three_param_baseline_bounded",
                      "PRODUCTION incumbent: three_param, b in [0.005,0.04]", fit_one)
