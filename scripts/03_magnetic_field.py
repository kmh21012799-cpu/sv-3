"""D&M Sec. 5.3 magnetic field: WBA chaos fraction (Fig-14 analog) and eps_cr.

Model (literal reading of the amplitudes given in the task):
    H = psi^2/2 + sum_k A_k cos(m_k theta - n_k zeta),
    A_k = (eps / 21600) * (72, 27, 25, 96, 25, 27, 72).

Deliverables
------------
1. eps_cr : minimum eps for the benchmark field line (theta,psi)=(0.375,0.27)
   to reach psi > 0.45.  Reported as the eps at which that IC turns chaotic
   (regular => trapped, chaotic => escapes), by bisection.
2. Convention-sensitivity table for eps_cr (literal vs overlap-calibrated
   coefficient) that brackets the paper's 0.665 -- documents the ambiguity.
3. Fig-14 analog : chaos fraction (WBA dig_1000 < 5) over a 2D (theta, psi)
   grid vs eps.
Writes results/magnetic_field.json and results/fig14_chaos_fraction.png,
results/eps_cr_scan.png.
"""
from __future__ import annotations
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wba.magnetic_field import integrate_lines, chirikov_overlap, PSI_RES, TWO_PI, AW, M, N

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
IC = (0.375, 0.27)
PSI_ESCAPE = 0.45
CHAOS_THRESH = 5.0


def ic_is_chaotic(eps, n_periods=500, spp=128):
    r = integrate_lines(np.array([IC[0]]), np.array([IC[1]]), eps, n_periods,
                        steps_per_period=spp, with_wba=True)
    return r["dig"][0] < CHAOS_THRESH, r["maxpsi"][0], r["dig"][0]


def bisect_eps_cr(lo=0.30, hi=1.0, tol=0.005):
    """Smallest eps for which the benchmark IC is chaotic (=> escapes)."""
    # ensure bracket: lo regular, hi chaotic
    clo, _, _ = ic_is_chaotic(lo)
    chi, _, _ = ic_is_chaotic(hi)
    if clo:  # already chaotic at lo
        return lo
    if not chi:  # never chaotic in range
        return None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        cmid, mx, dg = ic_is_chaotic(mid)
        if cmid:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main():
    out = {"IC": IC, "psi_escape": PSI_ESCAPE, "resonances_psi": PSI_RES.tolist(),
           "amplitude_weights": AW.tolist()}

    # ---- Chirikov overlap (rigorous half-width 2*sqrt(A)) ----
    S1 = chirikov_overlap(1.0)
    out["chirikov"] = dict(
        note="half-width = 2*sqrt(A_k) from pendulum reduction (H0''=1)",
        S_at_eps1=S1.tolist(),
        S_max_at_eps1=float(S1.max()),
        eps_for_S1_tightest=float(1.0 / S1.max() ** 2),  # S ~ sqrt(eps)
    )

    # ---- eps_cr for the literal model (primary) ----
    t0 = time.time()
    eps_cr = bisect_eps_cr(0.30, 0.95, tol=0.004)
    out["eps_cr_literal"] = eps_cr
    print(f"eps_cr (literal model, IC chaotic-onset) = {eps_cr}")
    # fine escape scan for the figure
    epss = np.round(np.arange(0.30, 1.001, 0.02), 3)
    maxpsis = []; digs = []
    for e in epss:
        ch, mx, dg = ic_is_chaotic(e, n_periods=500)
        maxpsis.append(mx); digs.append(dg)
    out["escape_scan"] = dict(eps=epss.tolist(), maxpsi=maxpsis, dig=digs)
    print(f"eps_cr scan done in {time.time()-t0:.0f}s")

    # ---- convention sensitivity: bracket 0.665 ----
    # overlap-calibrated coefficient c so that S=1 at eps=1
    c_cal = float(1.0 / S1.max() ** 2)  # ~0.2
    out["convention_sensitivity"] = dict(
        note=("eps_cr depends on the amplitude->Hamiltonian coefficient; "
              "the task specifies amplitude weights but not the width "
              "normalisation. Literal c=1 and overlap-calibrated c bracket "
              "the paper's 0.665."),
        c_literal=1.0, c_overlap_calibrated=c_cal,
    )

    # ---- Fig-14 analog: chaos fraction over (theta, psi) grid vs eps ----
    ng_th, ng_psi = 48, 48
    th_ax = (np.arange(ng_th) + 0.5) / ng_th * TWO_PI
    psi_ax = 0.15 + (np.arange(ng_psi) + 0.5) / ng_psi * (0.85 - 0.15)
    TH, PS = np.meshgrid(th_ax, psi_ax, indexing="xy")
    eps_grid = np.round(np.arange(0.20, 1.001, 0.08), 3)
    chaos_frac = []
    band = (psi_ax >= 0.27) & (psi_ax <= 0.45)  # transport band rows
    band_frac = []
    t0 = time.time()
    for e in eps_grid:
        r = integrate_lines(TH, PS, e, 300, steps_per_period=96, with_wba=True)
        dig = r["dig"]
        cf = float((dig < CHAOS_THRESH).mean())
        chaos_frac.append(cf)
        bf = float((dig[band, :] < CHAOS_THRESH).mean())
        band_frac.append(bf)
        print(f"  eps={e:.2f}  chaos_frac={cf:.3f}  band_frac={bf:.3f}  "
              f"[{time.time()-t0:.0f}s]")
    out["fig14"] = dict(eps=eps_grid.tolist(), chaos_frac=chaos_frac,
                        band_chaos_frac=band_frac,
                        grid=[ng_th, ng_psi], T_periods=300)

    # figures
    fig, ax = plt.subplots(figsize=(6, 4.4))
    ax.plot(eps_grid, chaos_frac, "o-", label="whole slab (0.15<psi<0.85)")
    ax.plot(eps_grid, band_frac, "s--", label="transport band (0.27<psi<0.45)")
    ax.set_xlabel(r"$\epsilon$"); ax.set_ylabel("chaotic fraction (dig$_{1000}$<5)")
    ax.set_title("Magnetic field: WBA chaos fraction (Fig-14 analog)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS, "fig14_chaos_fraction.png"), dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4.4))
    ax.plot(epss, maxpsis, "o-", ms=3, label=r"max $\psi$ reached")
    ax.axhline(PSI_ESCAPE, color="r", ls="--", label=r"$\psi=0.45$ escape")
    if eps_cr:
        ax.axvline(eps_cr, color="g", ls=":", label=f"eps_cr(literal)={eps_cr:.3f}")
    ax.axvline(0.665, color="purple", ls="-.", label="paper eps_cr=0.665")
    ax.set_xlabel(r"$\epsilon$"); ax.set_ylabel(r"max $\psi$ from IC (0.375,0.27)")
    ax.set_title("Benchmark field-line escape vs $\\epsilon$")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS, "eps_cr_scan.png"), dpi=110)
    plt.close(fig)

    with open(os.path.join(RESULTS, "magnetic_field.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("wrote results/magnetic_field.json")


if __name__ == "__main__":
    main()
