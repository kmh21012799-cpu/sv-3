"""Weighted Birkhoff Average (WBA) for maps, and the dig_T convergence test.

Reference: Duignan & Meiss, Physica D 449 (2023) 133749.

Discrete (map) WBA
------------------
For an orbit x_0, x_1, ... of a map, and observable h,

    WB_T(h)(x_0) = (1/S) * sum_{n=0}^{T-1} g(n/T) h(x_n),
    S            = sum_{n=0}^{T-1} g(n/T),
    g(s)         = exp(-w/[s(1-s)])      (shape only; C cancels via /S).

Because we divide by S = sum of weights, the continuous normalisation constant
C plays no role in the map case (answer to pre-code Q2).

Convergence / chaos test (dig_T)
--------------------------------
Two independent windows of length T are used (pre-code Q3):

    A = WB_T(h) over x_0 ... x_{T-1}          (start at x_0)
    B = WB_T(h) over x_T ... x_{2T-1}         (start at x_T = phi_T(x_0))

    absdig_T = -log10 |A - B|
    reldig_T = -log10 |(A - B) / mean(|A|,|B|)|     (skip if mean ~ 0)
    dig_T    = max(absdig_T, reldig_T)

Regular (Diophantine) orbits: super-convergence, dig_T >> 1 (~10-15 at T=1000).
Chaotic orbits: slow convergence, dig_T small (< 5 at T=1000).
"""

from __future__ import annotations

import numpy as np


def weights(T: int, w: float = 1.0) -> np.ndarray:
    """Normalised discrete weights g(n/T)/S, n = 0..T-1 (sum to 1)."""
    n = np.arange(T, dtype=float)
    s = n / T
    g = np.zeros(T)
    inside = (s > 0.0) & (s < 1.0)
    ss = s[inside]
    g[inside] = np.exp(-w / (ss * (1.0 - ss)))
    S = g.sum()
    return g / S


def wba_from_series(hvals: np.ndarray, w: float = 1.0) -> float:
    """WB_T(h)(x_0) from a length-T series h(x_0..x_{T-1})."""
    T = len(hvals)
    return float(np.dot(weights(T, w), hvals))


def dig_from_two_windows(hvals: np.ndarray, T: int, w: float = 1.0):
    """Compute (dig_T, absdig, reldig, A, B) from a length-2T series.

    hvals must contain h(x_0), ..., h(x_{2T-1}).
    """
    wts = weights(T, w)
    A = float(np.dot(wts, hvals[:T]))
    B = float(np.dot(wts, hvals[T:2 * T]))
    diff = abs(A - B)
    # absolute digits
    if diff == 0.0:
        absdig = 16.0  # at/below double precision floor
    else:
        absdig = -np.log10(diff)
    # relative digits (unusable when the mean is ~0)
    scale = 0.5 * (abs(A) + abs(B))
    if scale < 1e-12:
        reldig = -np.inf
    else:
        rel = diff / scale
        reldig = 16.0 if rel == 0.0 else -np.log10(rel)
    dig = max(absdig, reldig)
    # clip to a sane double-precision ceiling
    dig = min(dig, 16.0)
    return dig, absdig, reldig, A, B


def classify_orbit(step, x0, h, T: int, w: float = 1.0, chaos_thresh: float = 5.0):
    """Iterate `step` for 2T iterations, return dict with dig_T, rotation etc.

    Parameters
    ----------
    step : callable  x -> x_next   (single point, returns next point)
    x0   : initial condition (array-like)
    h    : callable  x -> float    (observable)
    T    : window length
    w    : weight parameter (fixed at 1.0 for the paper)
    chaos_thresh : dig_T below this => flagged chaotic.

    Returns
    -------
    dict with keys: dig, absdig, reldig, wba, is_regular, escaped
    """
    x = np.array(x0, dtype=float)
    hvals = np.empty(2 * T)
    escaped = False
    for n in range(2 * T):
        hvals[n] = h(x)
        if not np.all(np.isfinite(x)):
            escaped = True
            break
        x = step(x)
    if escaped:
        return dict(dig=0.0, absdig=0.0, reldig=-np.inf, wba=np.nan,
                    is_regular=False, escaped=True)
    dig, absdig, reldig, A, B = dig_from_two_windows(hvals, T, w)
    return dict(dig=dig, absdig=absdig, reldig=reldig, wba=A,
                is_regular=(dig >= chaos_thresh), escaped=False)
