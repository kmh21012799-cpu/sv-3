"""Chirikov-Taylor standard map, WBA sweeps, and Benettin Lyapunov exponents.

Map on the 2-torus (theta, p) in [0, 2pi)^2:

    p'     = p + K sin(theta)      (mod 2pi)
    theta' = theta + p'            (mod 2pi)

Reducing p mod 2pi does not change the torus dynamics because theta' shifts by
a multiple of 2pi.  On the torus, the K>2pi "accelerator" fixed points (which
accelerate p on the cylinder) appear as ordinary elliptic islands; WBA flags
them regular.

Observable
----------
Default classification observable is h = sin(theta): a smooth bounded function
on the torus.  For a Diophantine rotational circle *or* an island it is a
smooth quasi-periodic function of the orbit, so its Birkhoff average
super-converges; on chaotic orbits it does not.  (Using the raw action p as the
observable fails on accelerator islands, where p advances by 2pi per step on
the cylinder and the two windows differ by ~2pi*T; sin(theta) avoids this.)

The internal rotation number is obtained separately as the WBA of the wrapped
angle advance delta_n = wrap(theta_{n+1} - theta_n).
"""

from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


def sm_step(x, K: float):
    """Single-point standard-map step on the torus. x = [theta, p]."""
    theta, p = x
    p = (p + K * np.sin(theta)) % TWO_PI
    theta = (theta + p) % TWO_PI
    return np.array([theta, p])


def sm_step_batch(theta, p, K: float):
    """Vectorised standard-map step for arrays theta, p."""
    p = (p + K * np.sin(theta)) % TWO_PI
    theta = (theta + p) % TWO_PI
    return theta, p


def wrap_pi(x):
    """Wrap angle to (-pi, pi]."""
    return (x + np.pi) % TWO_PI - np.pi


def sweep_wba(theta0, p0, K: float, T: int, w: float = 1.0):
    """Vectorised WBA sweep over arrays of initial conditions.

    Returns dict of arrays (same shape as inputs):
      dig      : max(absdig, reldig) convergence digits
      wba_sin  : WB_T(sin theta) over window A
      rot      : WB_T(wrapped angle advance) = internal rotation number / step
    Uses observable h = sin(theta) for the dig test; window A = steps 0..T-1,
    window B = steps T..2T-1.
    """
    theta = np.array(theta0, dtype=float).ravel().copy()
    p = np.array(p0, dtype=float).ravel().copy()
    N = theta.size

    n = np.arange(T, dtype=float)
    s = n / T
    g = np.zeros(T)
    inside = (s > 0) & (s < 1)
    g[inside] = np.exp(-1.0 / (s[inside] * (1 - s[inside]))) if w == 1.0 \
        else np.exp(-w / (s[inside] * (1 - s[inside])))
    g /= g.sum()

    A_sin = np.zeros(N)
    B_sin = np.zeros(N)
    A_rot = np.zeros(N)  # rotation-number WBA (window A only, for the value)

    theta_prev = theta.copy()
    for k in range(2 * T):
        hsin = np.sin(theta)
        if k < T:
            A_sin += g[k] * hsin
        else:
            B_sin += g[k - T] * hsin
        theta_next, p = sm_step_batch(theta, p, K)
        if k < T:
            dtheta = wrap_pi(theta_next - theta)
            A_rot += g[k] * dtheta
        theta = theta_next

    diff = np.abs(A_sin - B_sin)
    with np.errstate(divide='ignore'):
        absdig = -np.log10(diff)
    absdig[diff == 0.0] = 16.0
    scale = 0.5 * (np.abs(A_sin) + np.abs(B_sin))
    rel = np.where(scale > 1e-12, diff / np.maximum(scale, 1e-300), np.inf)
    with np.errstate(divide='ignore'):
        reldig = -np.log10(rel)
    reldig[rel == 0.0] = 16.0
    reldig[~np.isfinite(rel)] = -np.inf
    dig = np.maximum(absdig, reldig)
    dig = np.minimum(dig, 16.0)

    escaped = ~np.isfinite(theta) | ~np.isfinite(p)
    dig[escaped] = 0.0
    return dict(dig=dig.reshape(np.shape(theta0)),
                wba_sin=A_sin.reshape(np.shape(theta0)),
                rot=(A_rot / TWO_PI).reshape(np.shape(theta0)))


def benettin_lyapunov(x0, K: float, N: int, renorm_every: int = 1):
    """Largest Lyapunov exponent by the Benettin tangent-vector method.

    Returns lambda (per iteration).  Tangent map of the standard map:
        d p'     = d p + K cos(theta) d theta
        d theta' = d theta + d p'
    """
    theta, p = float(x0[0]), float(x0[1])
    v = np.array([1.0, 0.0])
    v /= np.linalg.norm(v)
    lyap_sum = 0.0
    count = 0
    for i in range(N):
        c = K * np.cos(theta)
        # tangent step (same order as the map)
        dp = v[1] + c * v[0]
        dtheta = v[0] + dp
        v = np.array([dtheta, dp])
        # base step
        p = (p + K * np.sin(theta)) % TWO_PI
        theta = (theta + p) % TWO_PI
        if (i + 1) % renorm_every == 0:
            nrm = np.linalg.norm(v)
            lyap_sum += np.log(nrm)
            v /= nrm
            count += renorm_every
    return lyap_sum / count


def benettin_lyapunov_batch(theta0, p0, K: float, N: int):
    """Vectorised largest Lyapunov exponent over arrays of ICs."""
    theta = np.array(theta0, dtype=float).ravel().copy()
    p = np.array(p0, dtype=float).ravel().copy()
    M = theta.size
    vt = np.ones(M)
    vp = np.zeros(M)
    nrm = np.hypot(vt, vp)
    vt /= nrm
    vp /= nrm
    lyap = np.zeros(M)
    for i in range(N):
        c = K * np.cos(theta)
        dp = vp + c * vt
        dtheta = vt + dp
        vt, vp = dtheta, dp
        p = (p + K * np.sin(theta)) % TWO_PI
        theta = (theta + p) % TWO_PI
        nrm = np.hypot(vt, vp)
        lyap += np.log(nrm)
        vt /= nrm
        vp /= nrm
    return (lyap / N).reshape(np.shape(theta0))
