"""D&M (2023) Sec. 5.3 field-line Hamiltonian and WBA via an augmented ODE.

Field-line Hamiltonian (zeta is the time-like toroidal angle):

    H(theta, psi, zeta) = psi^2 / 2 + sum_k A_k cos(m_k theta - n_k zeta)

with the Farey-tree level-3 resonances and amplitudes

    (m,n) : (4,1) (3,1) (5,2) (2,1) (5,3) (3,2) (4,3)
    A_k   = (eps/21600) * (72, 27, 25, 96, 25, 27, 72)

Resonance (m,n) sits at psi = n/m.  Amplitudes are tuned (per D&M) so that all
neighbouring resonance pairs reach Chirikov overlap simultaneously at eps = 1.

Hamilton's equations (zeta = time):
    dtheta/dzeta =  dH/dpsi   =  psi
    dpsi/dzeta   = -dH/dtheta =  sum_k A_k m_k sin(m_k theta - n_k zeta)

WBA is integrated *with the same integrator* as an extra ODE (paper's method):
    dW/dzeta = g_norm(zeta / T) * h(state),     h = psi,     WB = W(T)/T,
with g_norm normalised so int_0^1 g = 1.  Two consecutive windows of length T
give the dig_T convergence estimate.
"""
from __future__ import annotations
import numpy as np

TWO_PI = 2.0 * np.pi

# resonances and base amplitude weights
M = np.array([4, 3, 5, 2, 5, 3, 4], dtype=float)
N = np.array([1, 1, 2, 1, 3, 2, 3], dtype=float)
AW = np.array([72, 27, 25, 96, 25, 27, 72], dtype=float)   # /21600 * eps
PSI_RES = N / M   # 0.25, 0.333, 0.4, 0.5, 0.6, 0.667, 0.75

# normalisation constant for g (w=1): from wba.weight (int g = 1)
_C_W1 = 142.25037577713


def amplitudes(eps: float) -> np.ndarray:
    return (eps / 21600.0) * AW


def chirikov_overlap(eps: float):
    """Return per-neighbour-pair overlap ratios S = (dw_k + dw_{k+1})/gap.

    Single-resonance island half-width in psi: dpsi_k = 2*sqrt(A_k) (from the
    pendulum reduction with H0 = psi^2/2, |d2H0/dpsi2| = 1).
    """
    A = amplitudes(eps)
    dw = 2.0 * np.sqrt(A)
    gaps = np.diff(PSI_RES)
    S = (dw[:-1] + dw[1:]) / gaps
    return S


def deriv(theta, psi, zeta, A):
    """Vectorised RHS of the field-line ODE for arrays theta, psi (scalar zeta)."""
    # dpsi = sum_k A_k m_k sin(m_k theta - n_k zeta)
    dpsi = np.zeros_like(theta)
    for k in range(len(M)):
        dpsi += A[k] * M[k] * np.sin(M[k] * theta - N[k] * zeta)
    dtheta = psi
    return dtheta, dpsi


def integrate_lines(theta0, psi0, eps, n_periods, steps_per_period=256,
                    with_wba=False, w=1.0, track_maxpsi=False):
    """RK4-integrate an ensemble of field lines.

    Returns dict.  If with_wba: integrates two consecutive length-T windows
    (T = n_periods toroidal periods) and returns dig_T etc. (h = psi).
    If track_maxpsi: returns running max/min of psi over a single window.
    """
    theta = np.array(theta0, dtype=float).ravel().copy()
    psi = np.array(psi0, dtype=float).ravel().copy()
    A = amplitudes(eps)
    dz = TWO_PI / steps_per_period

    if track_maxpsi:
        total_steps = n_periods * steps_per_period
        maxpsi = psi.copy(); minpsi = psi.copy()
        zeta = 0.0
        for _ in range(total_steps):
            theta, psi, zeta = _rk4(theta, psi, zeta, dz, A)
            np.maximum(maxpsi, psi, out=maxpsi)
            np.minimum(minpsi, psi, out=minpsi)
        return dict(maxpsi=maxpsi.reshape(np.shape(theta0)),
                    minpsi=minpsi.reshape(np.shape(theta0)),
                    theta=theta, psi=psi)

    if with_wba:
        T = n_periods * TWO_PI            # window length in zeta
        total_steps = 2 * n_periods * steps_per_period
        WA = np.zeros_like(psi); WB = np.zeros_like(psi)
        maxpsi = psi.copy(); minpsi = psi.copy()
        zeta = 0.0
        for i in range(total_steps):
            # weight argument s in (0,1) within the current window
            in_first = zeta < T
            s = (zeta / T) if in_first else ((zeta - T) / T)
            g = 0.0
            if 0.0 < s < 1.0:
                g = _C_W1 * np.exp(-w / (s * (1.0 - s)))
            # RK4 step for state; WBA accumulated with trapezoid via midpoint psi
            h0 = psi
            theta, psi, zeta = _rk4(theta, psi, zeta, dz, A)
            # accumulate dW = g * h * dz using average of endpoints
            h_avg = 0.5 * (h0 + psi)
            if in_first:
                WA += g * h_avg * dz
            else:
                WB += g * h_avg * dz
            np.maximum(maxpsi, psi, out=maxpsi)
            np.minimum(minpsi, psi, out=minpsi)
        A_wba = WA / T
        B_wba = WB / T
        diff = np.abs(A_wba - B_wba)
        with np.errstate(divide='ignore'):
            absdig = -np.log10(diff)
        absdig[diff == 0.0] = 16.0
        scale = 0.5 * (np.abs(A_wba) + np.abs(B_wba))
        rel = np.where(scale > 1e-12, diff / np.maximum(scale, 1e-300), np.inf)
        with np.errstate(divide='ignore'):
            reldig = -np.log10(rel)
        reldig[~np.isfinite(rel)] = -np.inf
        dig = np.maximum(absdig, reldig)
        dig = np.minimum(dig, 16.0)
        bad = ~np.isfinite(psi)
        dig[bad] = 0.0
        return dict(dig=dig.reshape(np.shape(theta0)),
                    wba=A_wba.reshape(np.shape(theta0)),
                    maxpsi=maxpsi.reshape(np.shape(theta0)),
                    minpsi=minpsi.reshape(np.shape(theta0)))

    raise ValueError("set with_wba or track_maxpsi")


def _rk4(theta, psi, zeta, dz, A):
    dt1, dp1 = deriv(theta, psi, zeta, A)
    dt2, dp2 = deriv(theta + 0.5 * dz * dt1, psi + 0.5 * dz * dp1, zeta + 0.5 * dz, A)
    dt3, dp3 = deriv(theta + 0.5 * dz * dt2, psi + 0.5 * dz * dp2, zeta + 0.5 * dz, A)
    dt4, dp4 = deriv(theta + dz * dt3, psi + dz * dp3, zeta + dz, A)
    theta = theta + dz / 6.0 * (dt1 + 2 * dt2 + 2 * dt3 + dt4)
    psi = psi + dz / 6.0 * (dp1 + 2 * dp2 + 2 * dp3 + dp4)
    theta = theta % TWO_PI
    zeta = zeta + dz
    return theta, psi, zeta
