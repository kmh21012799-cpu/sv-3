"""Standard-map WBA deliverables.

Produces, for each K in {0.5, 0.9716, 1.2, 2.0, 3.0, 6.9}:
  * a dig_T map over the (theta, p) torus  -> results/dig_map_K*.png
  * chaos fraction (dig_1000 < 5) on a 2D grid and on the theta0=0 line
  * WBA-vs-Lyapunov agreement + timing on a shared subsample

Writes results/standard_map_stats.json.
"""
from __future__ import annotations
import json, time, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wba.standard_map import sweep_wba, benettin_lyapunov_batch, TWO_PI

K_LIST = [0.5, 0.9716, 1.2, 2.0, 3.0, 6.9]
T = 1000
CHAOS_THRESH = 5.0
GRID = 360          # GRID x GRID torus grid for the dig maps
RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
os.makedirs(RESULTS, exist_ok=True)


def chunked_sweep(TH, PP, K, T, chunk=40000):
    """Run sweep_wba in flat chunks to bound memory; return dig, rot arrays."""
    th = TH.ravel(); pp = PP.ravel()
    dig = np.empty(th.size); rot = np.empty(th.size)
    for i in range(0, th.size, chunk):
        sl = slice(i, i + chunk)
        r = sweep_wba(th[sl], pp[sl], K, T)
        dig[sl] = r["dig"]; rot[sl] = r["rot"]
    return dig.reshape(TH.shape), rot.reshape(TH.shape)


def main():
    stats = {"T": T, "chaos_thresh": CHAOS_THRESH, "grid": GRID, "per_K": {}}
    # avoid the theta0=0 symmetry line and exact resonant p; small offset grid
    ax = (np.arange(GRID) + 0.5) / GRID * TWO_PI
    TH, PP = np.meshgrid(ax, ax, indexing="xy")

    for K in K_LIST:
        t0 = time.time()
        dig, rot = chunked_sweep(TH, PP, K, T)
        dt = time.time() - t0
        regular = dig >= CHAOS_THRESH
        chaos_frac_grid = float(1.0 - regular.mean())

        # theta0 = 0 line control (dominant symmetry line): p varies, theta=0
        pax = (np.arange(GRID) + 0.5) / GRID * TWO_PI
        digs_line, _ = chunked_sweep(np.zeros_like(pax), pax, K, T)
        chaos_frac_line = float((digs_line < CHAOS_THRESH).mean())

        # boundary counts (near-boundary orbits should be rare)
        near = np.sum((dig > 4.0) & (dig < 6.0))

        stats["per_K"][str(K)] = dict(
            chaos_frac_grid=chaos_frac_grid,
            chaos_frac_line=chaos_frac_line,
            regular_frac_grid=float(regular.mean()),
            n_boundary_4to6=int(near),
            n_total=int(dig.size),
            dig_median_regular=float(np.median(dig[regular])) if regular.any() else None,
            dig_median_chaotic=float(np.median(dig[~regular])) if (~regular).any() else None,
            sweep_seconds=dt,
        )
        print(f"K={K:6.4f}  chaos(grid)={chaos_frac_grid:.4f}  "
              f"chaos(line)={chaos_frac_line:.4f}  "
              f"bdry(4<dig<6)={near}  t={dt:.1f}s")

        # ---- dig map figure ----
        fig, ax2 = plt.subplots(figsize=(5.6, 5.0))
        im = ax2.imshow(np.clip(dig, 0, 16), origin="lower",
                        extent=[0, TWO_PI, 0, TWO_PI], aspect="auto",
                        cmap="viridis", vmin=0, vmax=16)
        ax2.set_title(f"Standard map  K={K}   WBA dig$_{{{T}}}$  (h=sin$\\theta$)")
        ax2.set_xlabel(r"$\theta$"); ax2.set_ylabel(r"$p$")
        cb = fig.colorbar(im, ax=ax2); cb.set_label(r"dig$_T$ (convergence digits)")
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS, f"dig_map_K{K}.png"), dpi=110)
        plt.close(fig)

    # ---- WBA vs Lyapunov agreement (subsample across all K at K=6.9 richest) ----
    Kc = 6.9
    sub = 120
    axs = (np.arange(sub) + 0.5) / sub * TWO_PI
    STH, SPP = np.meshgrid(axs, axs, indexing="xy")
    t0 = time.time()
    dig_s, _ = chunked_sweep(STH, SPP, Kc, T)
    t_wba = time.time() - t0
    t0 = time.time()
    lyap = benettin_lyapunov_batch(STH, SPP, Kc, 2 * T)
    t_lyap = time.time() - t0
    wba_reg = (dig_s >= CHAOS_THRESH).ravel()
    lyap_reg = (lyap.ravel() < 0.03)   # near-zero Lyapunov => regular
    agree = float(np.mean(wba_reg == lyap_reg))
    stats["wba_vs_lyapunov"] = dict(
        K=Kc, n=int(wba_reg.size), agreement=agree,
        wba_regular_frac=float(wba_reg.mean()),
        lyap_regular_frac=float(lyap_reg.mean()),
        t_wba_seconds=t_wba, t_lyap_seconds=t_lyap,
        speedup_lyap_over_wba=t_lyap / t_wba,
    )
    print(f"\nWBA vs Lyapunov @K={Kc}: agreement={agree:.4f}  "
          f"t_wba={t_wba:.1f}s t_lyap={t_lyap:.1f}s")

    # scatter figure
    fig, axf = plt.subplots(figsize=(5.4, 4.6))
    sc = axf.scatter(lyap.ravel(), dig_s.ravel(), s=4, alpha=0.35, c=wba_reg,
                     cmap="coolwarm")
    axf.axhline(CHAOS_THRESH, color="k", ls="--", lw=0.8, label="dig=5")
    axf.axvline(0.03, color="gray", ls=":", lw=0.8, label="lyap=0.03")
    axf.set_xlabel("Benettin Lyapunov exponent"); axf.set_ylabel("WBA dig$_{1000}$")
    axf.set_title(f"WBA vs Lyapunov, standard map K={Kc}")
    axf.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS, "wba_vs_lyapunov.png"), dpi=110)
    plt.close(fig)

    with open(os.path.join(RESULTS, "standard_map_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    print("wrote results/standard_map_stats.json")


if __name__ == "__main__":
    main()
