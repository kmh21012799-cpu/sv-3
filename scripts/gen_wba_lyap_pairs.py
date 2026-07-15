"""Regenerate and SAVE the per-point (dig, lambda) pairs for the WBA-Lyapunov
cross-check at K=6.9.

These 14,400 pairs were previously rendered straight to PNG and not stored, so
restyling the figure forced a recompute. This script saves them to
results/wba_lyap_pairs.npz so that never happens again.

The computation is identical to the WBA-vs-Lyapunov block in
scripts/01_standard_map.py (K=6.9, 120x120 grid, T=1000, dig via sweep_wba,
Benettin largest Lyapunov via 2T iterations). It is deterministic (no RNG), so
it reproduces the summary already recorded in results/standard_map_stats.json.
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from wba.standard_map import sweep_wba, benettin_lyapunov_batch, TWO_PI

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")

Kc = 6.9
T = 1000
SUB = 120                # 120 x 120 = 14400 initial conditions
CHAOS_THRESH = 5.0       # dig >= 5  => regular (WBA)
LYAP_THRESH = 0.03       # lambda < 0.03 => regular (Benettin)


def main():
    axs = (np.arange(SUB) + 0.5) / SUB * TWO_PI
    STH, SPP = np.meshgrid(axs, axs, indexing="xy")

    # dig from the weighted Birkhoff sweep (h = sin theta), same as script 01
    r = sweep_wba(STH.ravel(), SPP.ravel(), Kc, T)
    dig = r["dig"].astype(float)

    # largest Lyapunov exponent, Benettin, 2T iterations
    lam = benettin_lyapunov_batch(STH, SPP, Kc, 2 * T).ravel().astype(float)

    theta0 = STH.ravel().astype(float)
    p0 = SPP.ravel().astype(float)

    wba_reg = dig >= CHAOS_THRESH
    lyap_reg = lam < LYAP_THRESH
    agreement = float(np.mean(wba_reg == lyap_reg))

    out = os.path.join(RESULTS, "wba_lyap_pairs.npz")
    np.savez_compressed(
        out,
        dig=dig, lam=lam, theta0=theta0, p0=p0,
        K=np.float64(Kc), T=np.int64(T), sub=np.int64(SUB),
        chaos_thresh=np.float64(CHAOS_THRESH), lyap_thresh=np.float64(LYAP_THRESH),
    )

    # cross-check against the recorded summary
    both_reg = int(np.sum(wba_reg & lyap_reg))
    both_chaos = int(np.sum(~wba_reg & ~lyap_reg))
    disagree = int(np.sum(wba_reg != lyap_reg))
    print(f"saved {out}")
    print(f"n = {dig.size}")
    print(f"agreement          = {agreement:.5f}  (recorded 0.99986)")
    print(f"WBA regular frac   = {wba_reg.mean():.4f}  (recorded 0.0043)")
    print(f"Lyap regular frac  = {lyap_reg.mean():.4f}  (recorded 0.0042)")
    print(f"both regular       = {both_reg}")
    print(f"both chaotic       = {both_chaos}")
    print(f"disagree           = {disagree}")
    print(f"dig  range         = [{dig.min():.2f}, {dig.max():.2f}]")
    print(f"lam  range         = [{lam.min():.4f}, {lam.max():.4f}]")


if __name__ == "__main__":
    main()
