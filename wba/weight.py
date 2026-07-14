"""Weighted-Birkhoff bump weight function.

Duignan & Meiss (2023, Physica D 449, 133749), Eq. defining the weight:

    g_w(s) = C * exp( -w / [ s (1-s) ] ),   s in (0,1)
           = 0                              otherwise

The constant C normalises the *continuous* integral to one:

    C = 1 / integral_0^1 exp(-w/[s(1-s)]) ds .

For w = 1 the paper quotes C ~= 142.2503758.

Notes
-----
* For the *continuous* (flow) case the WBA is computed as an extra ODE
  dW/dt = g_w(t/T) h(phi_t x0); there the normalisation constant matters and
  we use the properly normalised g_w so that WB_T = W(T)/T with int g = 1.
* For the *map* (discrete) case the WBA is a normalised weighted sum and the
  constant C cancels (we divide by the sum of the weights), so only the
  *shape* exp(-w/[s(1-s)]) is needed.  See wba.birkhoff.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate


def bump_shape(s, w: float = 1.0):
    """Unnormalised bump exp(-w/[s(1-s)]) on (0,1), zero elsewhere.

    Vectorised and safe at the endpoints (returns 0 there).
    """
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    inside = (s > 0.0) & (s < 1.0)
    ss = s[inside]
    out[inside] = np.exp(-w / (ss * (1.0 - ss)))
    return out


def normalisation_constant(w: float = 1.0) -> float:
    """C = 1 / int_0^1 exp(-w/[s(1-s)]) ds, by high-accuracy quadrature."""
    val, err = integrate.quad(lambda s: np.exp(-w / (s * (1.0 - s))),
                              0.0, 1.0, limit=400, epsabs=1e-14, epsrel=1e-14)
    return 1.0 / val


def g_w(s, w: float = 1.0, C: float | None = None):
    """Normalised weight g_w(s) with int_0^1 g_w = 1."""
    if C is None:
        C = normalisation_constant(w)
    return C * bump_shape(s, w)


# Cache the w=1 constant used throughout the paper.
C_W1 = normalisation_constant(1.0)


if __name__ == "__main__":
    C = normalisation_constant(1.0)
    print(f"C(w=1)                = {C:.10f}")
    print(f"paper value           = 142.2503758")
    print(f"abs difference        = {abs(C - 142.2503758):.3e}")
    # Verify int g = 1 with the normalised weight.
    val, _ = integrate.quad(lambda s: g_w(s, 1.0, C), 0.0, 1.0,
                            limit=400, epsabs=1e-14, epsrel=1e-14)
    print(f"int_0^1 g_1(s) ds     = {val:.12f}")
    for w in (0.5, 2.0, 3.0):
        print(f"C(w={w})              = {normalisation_constant(w):.10f}")
