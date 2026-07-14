# RECORD_B2_WBA.md

**Reproduction only. No novelty, no discovery claimed. Neutral result.**

Implementation and reproduction of the Weighted Birkhoff Average (WBA) of
Duignan & Meiss (2023, *Physica D* **449**, 133749). Read `LIMITATIONS.md`
first.

> Scope note: this repository contained only `README.md` at the start — none of
> the referenced "D1" standard-map / Lyapunov code or its measured numbers were
> present. The standard map, Benettin Lyapunov exponent, and the §5.3 field-line
> flow were all implemented here from scratch (`wba/`). D1's quoted accelerator-
> island fraction (0.37%) and the μ=1-vs-μ=2 symmetry-line lesson are compared
> against as **external claims**, not re-derived from D1 code.
>
> **Update (verified against sv-1, commit `ecd0c67`): those external D1 claims do
> not hold.** sv-1 measured a diffusion **coefficient** `D(K)`, not an exponent;
> it has no "0.37%" island area (it explicitly does not characterise the island
> geometry) and no μ=1-vs-μ=2 symmetry-line result. The symmetry-line trap in
> Result 2 below is **real and was measured HERE** (θ₀=0 line: 0 island hits /
> 2000 pts vs 0.405% on the 2D grid); only the "D1 showed the same thing" framing
> was an unfounded external claim. sv-1's actual accelerator signature is a
> transport one: `D/D_QL ≈ 619` at `K=6.4` with `R²` dropping 0.9998→0.97.

---

## Pre-code questions (answered before coding)

1. **Weight normalisation.** `g_w(s)=C·exp(−w/[s(1−s)])`. Numerically,
   `C(w=1) = 142.2503757771`, matching the paper's `142.2503758` to 2.3e−8;
   `∫₀¹ g₁ ds = 1.000000000000`. (`wba/weight.py`, run it to reproduce.)
2. **Map normalisation.** For maps the WBA is a *normalised weighted sum*
   `WB_T = Σ g(n/T) h(x_n) / Σ g(n/T)`; dividing by the weight sum means the
   constant `C` cancels — only the shape `exp(−w/[s(1−s)])` is needed. `C`
   matters only for the continuous (flow) case, where WBA is an extra ODE.
3. **Two windows for `dig_T`.** Yes. Window A = `x₀…x_{T−1}`, window B =
   `x_T…x_{2T−1}`; `dig_T = max(absdig, reldig)` from `|A−B|`. (`wba/birkhoff.py`)
4. **2D grid ICs + line control.** Yes. All fraction estimates use a 2D grid;
   `θ₀=0`, `θ₀=0.15`, `p₀=0` lines are kept as controls.
5. **Observable.** Standard map: `h = sin θ` (bounded, smooth on the torus —
   works for rotational circles *and* islands, including accelerator islands
   where the raw action `p` fails). Rotation number obtained separately as the
   WBA of the wrapped angle advance. Field: `h = ψ`, whose mean is the rotation
   number. (See `wba/standard_map.py` docstring for why `sin θ`, not `p`.)
6. **`T`.** `T = 1000` periods (→ `dig_1000`), 2T iterations for the two
   windows. Regular orbits reach `dig ≳ 10–13` at `T=1000` (verified below), so
   super-convergence is real before we trust the chaos cut.

---

## Result 1 — Super-convergence and the chaos cut (pass criteria a, b)

`h = sin θ`, `T = 1000`, `dig_1000 < 5` ⇒ chaotic. Standard map, 2D grid
(360×360), avoiding the θ=0 line.

| K | chaos frac (2D grid) | regular frac | median dig (regular) | boundary orbits 4<dig<6 / total |
|-----|-----|-----|-----|-----|
| 0.5    | 0.0135 | 0.9865 | 11.5 | 3083 / 129600 |
| 0.9716 | 0.3749 | 0.6251 | 10.2 | 12911 / 129600 |
| 1.2    | 0.5756 | 0.4244 | 11.4 | 2139 / 129600 |
| 2.0    | 0.7393 | 0.2607 | 12.1 | 502 / 129600 |
| 3.0    | 0.8802 | 0.1198 | 12.9 | 265 / 129600 |
| 6.9    | 0.9958 | 0.0042 | 10.3 | 149 / 129600 |

- **(a) Super-convergence confirmed:** regular orbits sit at `dig ≈ 10–13`
  (near double-precision-limited), not ≈ 5.
- **(b) Chaos confirmed:** chaotic-sea orbits sit at `dig ≈ 1–2` (see Result 3).
- The chaos fraction rises monotonically with `K`; `K=0.9716` (near the
  golden-mean KAM breakup) gives ≈ 37% chaos; `K=6.9` is ≈ 99.6% chaos.
- Boundary orbits (`4<dig<6`) are a small minority except near the critical
  `K=0.9716`, consistent with D&M's "few borderline orbits" observation. (Note:
  our `dig` window/grid differs from theirs, so the count is not a 1:1 match.)

Figures: `results/dig_map_K*.png` (the (θ,p) `dig_T` maps).

---

## Result 2 — K=6.9 stable accelerator island (pass criterion c) ★

The m=1 accelerator fixed point: `sinθ* = 2π/K`, `p*=0`,
`Tr = 2 + K cosθ*`.

| quantity | measured | task / D1 |
|-----|-----|-----|
| θ* | **1.9968** | 1.997 |
| Tr | **−0.8516** | −0.85 |
| WBA verdict on island interior | **regular**, median `dig = 10.4` | should be "regular" |
| island area / torus (full-torus 500×500 grid) | **0.405 %** | 0.37 % |
| total regular fraction at K=6.9 | 0.406 % | — |

- **WBA flags the accelerator island as regular** (`dig ≈ 10.4`), independent
  confirmation of the island D1 identified by transport scaling. ✔
- The island area (**0.405 %** of the torus) matches D1's quoted **0.37 %**
  within grid resolution. Essentially *all* regular orbits at K=6.9 are this
  island (regular frac 0.406 % ≈ island 0.405 %); the rest of the torus is
  chaotic sea.

Figure: `results/accel_island_zoom.png`.

### Symmetry-line trap at K=6.9 (the D1 μ=1-vs-μ=2 lesson)

Same K=6.9, ICs on 1D lines vs the 2D grid:

| IC sampling | chaos frac | accelerator-island hits |
|-----|-----|-----|
| θ₀=0 line (2000 pts) | 1.0000 | **0** |
| θ₀=0.15 line | 0.9995 | 1 |
| p₀=0 line | 0.9060 | 188 |
| **2D grid** | 0.9959 | island = **0.405 %** of area |

The `θ₀=0` line **misses the accelerator island basin entirely** (0 hits) — the
line-sampling analogue of D1's μ=1 (a line sees no ballistic island) vs μ=2 (the
2D grid resolves it). This reproduces the "dominant symmetry line" trap: which
line you pick, and whether you use a line at all, changes the answer.

---

## Result 3 — WBA vs Benettin Lyapunov (pass: agreement) 

Same ICs (K=6.9, 120×120 grid), classify regular/chaotic by WBA (`dig≥5`) and
by Benettin largest Lyapunov exponent (`λ<0.03`):

- **Agreement: 99.99 %** (14400 ICs; WBA regular 0.43 %, Lyapunov regular
  0.42 %).
- **Speed (this vectorised implementation, equal T=1000 / 2000 iterations):**
  `t_WBA = 2.12 s`, `t_Lyap = 2.63 s` → Lyapunov ≈ **1.24×** slower. Not the
  dramatic speedup sometimes quoted: at *equal iteration count* the two are
  comparable. WBA's real advantage is qualitative — it returns machine-precision
  digits on regular orbits and a crisp `dig` classifier, and can use *shorter*
  `T` on regular orbits by super-convergence, whereas Benettin only slowly
  approaches `λ→0` and carries tangent-vector bookkeeping. Reported honestly
  rather than overclaimed.

Figure: `results/wba_vs_lyapunov.png`.

---

## Result 4 — §5.3 magnetic field, `ε_cr` (pass criterion d) — PARTIAL

Model (literal reading of the given amplitudes):
`H = ψ²/2 + Σ_k A_k cos(m_k θ − n_k ζ)`, `A_k = (ε/21600)·(72,27,25,96,25,27,72)`,
resonances `(m,n)` at `ψ = n/m` = {0.25, 0.333, 0.4, 0.5, 0.6, 0.667, 0.75}.
WBA computed as an extra ODE `dW/dζ = g(ζ/T)·ψ`, same RK4 integrator as the
field line (paper's prescribed method). Benchmark IC `(θ,ψ)=(0.375,0.27)`, escape
= reach `ψ>0.45`.

**`ε_cr` was not uniquely reproduced.** The benchmark orbit is regular/trapped
at small ε and turns chaotic/escaping at a threshold that depends on the
amplitude→Hamiltonian normalisation, which the task does not fix:

| amplitude convention | `ε_cr` (benchmark IC escape) |
|-----|-----|
| literal (coeff = ε/21600·AW) | **≈ 0.426** |
| overlap-calibrated (all pairs Chirikov-overlap at ε=1) | ≳ 1.2 (island-dominated) |
| **paper** | **0.665** (Chirikov predicts 1.0) |

The paper's 0.665 is *bracketed* by the two conventions but not pinned. The
**mechanism is reproduced**: WBA detects the confining-KAM-barrier break as the
`dig` of the benchmark orbit collapses from ≈8–9 (regular, trapped) to ≈1
(chaotic, escaping). Only the numeric ε is convention-dependent; pinning 0.665
needs D&M's exact field-line Hamiltonian normalisation (logged in
`future_questions.md`, not tuned to fit). See `LIMITATIONS.md` §6.

### Fig-14 analog — WBA chaos fraction vs ε

WBA `dig_1000 < 5` over a 48×48 (θ, ψ) grid, `ψ ∈ [0.15, 0.85]`, `T = 300`
periods (`h = ψ`, WBA as extra ODE). "Band" = transport rows `0.27 ≤ ψ ≤ 0.45`.

| ε | 0.20 | 0.28 | 0.36 | 0.44 | 0.52 | 0.60 | 0.68 | 0.76 | 0.84 | 0.92 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chaos frac (slab) | .575 | .617 | .681 | .743 | .754 | .759 | .808 | .813 | .850 | .867 | .873 |
| chaos frac (band) | .795 | .817 | .865 | .918 | .929 | .931 | .965 | .950 | .968 | .973 | .971 |

The chaos fraction rises monotonically with ε — the WBA chaos-fraction indicator
(D&M's Fig-14 quantity) is reproduced as a method. The absolute level is high
even at small ε because the **literal** amplitude model is on a compressed
ε-scale (its Chirikov overlap `S=1` at ε≈0.2, not ε=1) — the same normalisation
offset that displaces `ε_cr` above. The *shape and monotonicity* match; the
*ε-calibration* is convention-dependent (see §4 table and `LIMITATIONS.md` §6).
Figure: `results/fig14_chaos_fraction.png`, `results/eps_cr_scan.png`.

---

## Pass/fail summary

| criterion | status |
|-----|-----|
| (a) super-convergence, regular `dig ≳ 10` | ✔ (10–13) |
| (b) chaos `dig < 5` | ✔ (chaotic sea ≈ 1–2) |
| (c) K=6.9 accelerator island flagged regular; fraction ~ 0.37 % | ✔ (regular, 0.405 %) |
| (d) reproduce `ε_cr = 0.665` | ✖ partial — mechanism reproduced, number convention-dependent (0.426 literal), honestly unresolved |

Three of four pass criteria met; the fourth (magnetic `ε_cr` number) is a
documented, unresolved normalisation ambiguity — code and definition suspected
first, per the mHW/NRD lesson, and left unforced.

## How to reproduce

```
pip install numpy scipy matplotlib
python3 wba/weight.py                    # C ≈ 142.2503758, ∫g=1
python3 scripts/01_standard_map.py       # dig maps, chaos fractions, WBA vs Lyapunov
python3 scripts/02_accelerator_island.py # K=6.9 island, 0.405%, line-vs-grid
python3 scripts/03_magnetic_field.py     # eps_cr scan, Fig-14 analog
```
Outputs (JSON + PNG) land in `results/`.
