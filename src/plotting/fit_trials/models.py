"""Model functions shared across JER fit trials.

Convention: E in GeV, all models return sigma/E. Parameter names use physics
letters (a=stochastic ~1/sqrt(E), b=constant, c=noise ~1/E, plus extras) - NOT
the swapped A/B/C display labels used elsewhere in the dashboard.
"""
import numpy as np


def m_stoch_const_noise(E, a, b, c):
    """Classic linear sum: a/sqrt(E) + b + c/E."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E


def m_stoch_const(E, a, b):
    """a/sqrt(E) + b."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b


def m_quad3(E, a, b, c):
    """Quadrature sum: sqrt((a/sqrt(E))^2 + b^2 + (c/E)^2)."""
    E = np.asarray(E, float)
    return np.sqrt((a / np.sqrt(E)) ** 2 + b ** 2 + (c / E) ** 2)


def m_quad2(E, a, b):
    """Quadrature sum: sqrt((a/sqrt(E))^2 + b^2)."""
    E = np.asarray(E, float)
    return np.sqrt((a / np.sqrt(E)) ** 2 + b ** 2)


def m_logE(E, a, b, c, d):
    """a/sqrt(E) + b + c/E + d*log(E)  (log term for slow high-E drift)."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E + d * np.log(E)


def m_stoch_const_logE(E, a, b, d):
    """a/sqrt(E) + b + d*log(E)  (drop noise, add log)."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + d * np.log(E)


def m_inv_sqrt_log(E, a, b, d):
    """a/sqrt(E) + b + d*log(E)/E  (log-suppressed noise)."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + d * np.log(E) / E


def m_power_const(E, a, p, b):
    """a * E^(-p) + b  (free exponent generalising the 1/sqrt(E) term)."""
    E = np.asarray(E, float)
    return a * E ** (-p) + b


def m_power_noise(E, a, b, c, p):
    """a/sqrt(E) + b + c * E^(-p)  (free exponent on the noise-like term)."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c * E ** (-p)


def m_stoch_const_noise2(E, a, b, c, e):
    """a/sqrt(E) + b + c/E + e/E^2  (extra steep low-E term)."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E + e / E ** 2


def m_two_power(E, a, p, c, q):
    """a*E^(-p) + c*E^(-q)  (sum of two power laws, no constant)."""
    E = np.asarray(E, float)
    return a * E ** (-p) + c * E ** (-q)


# --- Round 2 models (informed by round 1: e/E^2 and log(E) were the winners) ---

def m_logE_invE2(E, a, b, c, e):
    """a/sqrt(E) + b + c/E + e/E^2, but seeded/used as the log+steep combo's peer.

    NOTE: identical algebra to m_stoch_const_noise2; kept as a named alias so
    round-2 trials read clearly. Round 1 trial_20 established this as the leader.
    """
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E + e / E ** 2


def m_logE_and_invE2(E, a, b, c, e, d):
    """a/sqrt(E) + b + c/E + e/E^2 + d*log(E)  (both round-1 winners together). 5 params."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E + e / E ** 2 + d * np.log(E)


def m_invE15(E, a, b, c):
    """a/sqrt(E) + b + c/E^1.5  (steeper-than-noise low-E term, still 3 params)."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E ** 1.5


def m_invE15_noise(E, a, b, c, e):
    """a/sqrt(E) + b + c/E + e/E^1.5  (noise + a slightly steeper term). 4 params."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E + e / E ** 1.5


def m_free_low_power(E, a, b, c, r):
    """a/sqrt(E) + b + c/E^r  (free steep exponent r, nominally in [1,3]). 4 params."""
    E = np.asarray(E, float)
    return a / np.sqrt(E) + b + c / E ** r


def m_power_const_noise(E, a, p, b, c):
    """a*E^(-p) + b + c/E  (free stochastic exponent + constant + noise). 4 params."""
    E = np.asarray(E, float)
    return a * E ** (-p) + b + c / E


def m_logE_invE2_noStoch(E, b, c, e, d):
    """b + c/E + e/E^2 + d*log(E)  (drop the 1/sqrt(E) stochastic term). 4 params.

    Tests whether the sqrt term is even needed once steep + log terms are present.
    """
    E = np.asarray(E, float)
    return b + c / E + e / E ** 2 + d * np.log(E)
