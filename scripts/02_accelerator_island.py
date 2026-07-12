"""K=6.9 stable accelerator island: does WBA flag it as regular?

The m=1 accelerator fixed point of the standard map satisfies
    sin(theta*) = 2 pi / K,   p* = 0   (mod 2 pi).
For K = 6.9:  theta* = pi - arcsin(2 pi/K) ~= 2.000,  and the trace
    Tr = 2 + K cos(theta*) ~= -0.87   (|Tr| < 2  =>  elliptic / stable).
On the cylinder this fixed point accelerates p by 2 pi per step (ballistic
transport); on the torus it is an ordinary elliptic island.  WBA should flag
the island interior as *regular* (super-convergent), which is the core check:
an independent confirmation of the island D1 found by transport scaling.

Deliverables:
  * WBA dig on a fine zoom around the fixed point            -> results/accel_island_zoom.png
  * island area fraction (of the 4 pi^2 torus) vs the quoted 0.37%
  * a 1D IC line vs 2D grid demonstration that a line can miss the island
Writes results/accelerator_island.json
"""
from __future__ import annotations
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wba.standard_map import sweep_wba, TWO_PI

K = 6.9
T = 1000
CHAOS_THRESH = 5.0
RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


def chunked(th, pp, K, T, chunk=40000):
    th = th.ravel(); pp = pp.ravel()
    dig = np.empty(th.size); rot = np.empty(th.size); wsin = np.empty(th.size)
    for i in range(0, th.size, chunk):
        sl = slice(i, i + chunk)
        r = sweep_wba(th[sl], pp[sl], K, T)
        dig[sl] = r["dig"]; rot[sl] = r["rot"]; wsin[sl] = r["wba_sin"]
    return dig, rot, wsin


def main():
    theta_star = np.pi - np.arcsin(TWO_PI / K)
    Tr = 2 + K * np.cos(theta_star)
    print(f"accelerator fixed point: theta*={theta_star:.4f}  p*=0  Tr={Tr:.4f}")

    out = dict(K=K, T=T, theta_star=float(theta_star), trace=float(Tr))

    # ---- 1) fine zoom around the fixed point to see + measure the island ----
    half = 0.9  # zoom half-width in theta and p about (theta*, 0)
    nz = 700
    tax = theta_star + np.linspace(-half, half, nz)
    pax = 0.0 + np.linspace(-half, half, nz)
    TH, PP = np.meshgrid(tax % TWO_PI, pax % TWO_PI, indexing="xy")
    dig, rot, wsin = chunked(TH, PP, K, T)
    dig = dig.reshape(nz, nz); rot = rot.reshape(nz, nz)

    regular = dig >= CHAOS_THRESH
    # Island membership: regular AND rotation number matches the island's.
    # The island's internal rotation clusters tightly; identify it as the
    # dominant regular cluster located around the fixed point.
    reg_rot = rot[regular]
    # island rot ~ rotation at points nearest the fixed point
    center_mask = ((TH.reshape(nz, nz) - theta_star) ** 2 +
                   (PP.reshape(nz, nz)) ** 2) < 0.02
    if (regular & center_mask).any():
        rot_island = np.median(rot[regular & center_mask])
    else:
        rot_island = np.median(reg_rot) if reg_rot.size else 0.0
    island = regular & (np.abs(rot - rot_island) < 0.02)

    cell_area = (2 * half / (nz - 1)) ** 2
    island_area = float(island.sum() * cell_area)
    island_frac = island_area / (TWO_PI ** 2)
    out["island"] = dict(
        n_regular_in_zoom=int(regular.sum()),
        n_island=int(island.sum()),
        rot_island=float(rot_island),
        island_area=island_area,
        island_area_fraction_of_torus=island_frac,
        island_area_fraction_pct=100 * island_frac,
        quoted_D1_pct=0.37,
        dig_median_island=float(np.median(dig[island])) if island.any() else None,
    )
    print(f"island: n={int(island.sum())} regular, area/torus = "
          f"{100*island_frac:.3f}%  (D1 quote 0.37%)  "
          f"median dig in island = {np.median(dig[island]):.1f}")

    # figure: dig zoom with island contour
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.6))
    im0 = ax[0].imshow(np.clip(dig, 0, 16), origin="lower",
                       extent=[tax[0], tax[-1], pax[0], pax[-1]],
                       aspect="auto", cmap="viridis", vmin=0, vmax=16)
    ax[0].plot(theta_star, 0, "rx", ms=10, label="accel. fixed point")
    ax[0].set_title(f"K={K}  WBA dig$_{{{T}}}$ (island = high dig)")
    ax[0].set_xlabel(r"$\theta$"); ax[0].set_ylabel(r"$p$"); ax[0].legend(fontsize=8)
    fig.colorbar(im0, ax=ax[0], label="dig$_T$")
    ax[1].imshow(island, origin="lower",
                 extent=[tax[0], tax[-1], pax[0], pax[-1]], aspect="auto",
                 cmap="Greys")
    ax[1].plot(theta_star, 0, "rx", ms=10)
    ax[1].set_title("identified accelerator island (regular + matching rot)")
    ax[1].set_xlabel(r"$\theta$"); ax[1].set_ylabel(r"$p$")
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS, "accel_island_zoom.png"), dpi=110)
    plt.close(fig)

    # ---- 2) island area fraction over the FULL torus (independent estimate) ----
    ng = 500
    ax_full = (np.arange(ng) + 0.5) / ng * TWO_PI
    THf, PPf = np.meshgrid(ax_full, ax_full, indexing="xy")
    digf, rotf, _ = chunked(THf, PPf, K, T)
    regf = digf >= CHAOS_THRESH
    islandf = regf & (np.abs(rotf - rot_island) < 0.02)
    frac_full = float(islandf.mean())
    out["island"]["island_frac_full_torus_pct"] = 100 * frac_full
    out["island"]["regular_frac_full_torus_pct"] = 100 * float(regf.mean())
    print(f"full-torus grid ({ng}x{ng}): island fraction = {100*frac_full:.3f}%  "
          f"total regular = {100*regf.mean():.3f}%")

    # ---- 3) line vs grid: can a 1D line miss the island basin? ----
    npts = 2000
    line_results = {}
    for name, (thL, pL) in {
        "p0=0_line": ((np.arange(npts) + 0.5) / npts * TWO_PI, np.zeros(npts)),
        "theta0=0_line": (np.zeros(npts), (np.arange(npts) + 0.5) / npts * TWO_PI),
        "theta0=0.15_line": (np.full(npts, 0.15), (np.arange(npts) + 0.5) / npts * TWO_PI),
    }.items():
        dL, rL, _ = chunked(np.asarray(thL, float), np.asarray(pL, float), K, T)
        regL = dL >= CHAOS_THRESH
        islL = regL & (np.abs(rL - rot_island) < 0.02)
        line_results[name] = dict(
            chaos_frac=float((dL < CHAOS_THRESH).mean()),
            n_island_hits=int(islL.sum()),
        )
        print(f"line {name:18s}: chaos_frac={line_results[name]['chaos_frac']:.4f}  "
              f"island hits={line_results[name]['n_island_hits']}")
    out["line_vs_grid"] = dict(
        grid_island_frac_pct=100 * frac_full,
        grid_chaos_frac=float(1 - regf.mean()),
        lines=line_results,
    )

    with open(os.path.join(RESULTS, "accelerator_island.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("wrote results/accelerator_island.json")


if __name__ == "__main__":
    main()
