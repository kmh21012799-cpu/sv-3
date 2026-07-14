# WBA — Weighted Birkhoff Average

**Dynamics axis of a summer toolkit.**
**Full narrative: https://github.com/kmh21012799-cpu/sv1_comp**

## What this repository does

Implements the Weighted Birkhoff Average (WBA) of Duignan & Meiss
(2023, *Physica D* **449**, 133749) as a classifier for regular versus chaotic
orbits, and reproduces its standard-map and magnetic-field benchmarks.

WBA measures how fast a time-average converges along an orbit. A regular
(Diophantine) orbit makes the weighted average super-converge; a chaotic orbit
does not. The convergence is read off as a digit count,

```
dig = max(absdig, reldig),   absdig = -log10 | WBA[0,T] - WBA[T,2T] |,
```

the number of digits to which two consecutive time windows agree. Here
`T = 1000` map periods and the two windows are `x_0..x_{T-1}` and
`x_T..x_{2T-1}`. Empirically, regular orbits reach `dig ≈ 10–13` and chaotic
orbits sit at `dig ≈ 1–2`; the classifier cut used throughout is

```
dig_1000 < 5   =>   chaotic.
```

The weight is the smooth bump

```
g_w(s) = C * exp( -w / [ s (1-s) ] ),   s in (0,1),
```

with `w = 1` fixed (no tuning). Its normalisation constant was checked against
the paper: `C(w=1) = 142.2503757771`, matching the reported `142.2503758` to
2.3e-8, with `∫_0^1 g_1 ds = 1`. For maps the WBA is a normalised weighted sum,
so `C` cancels; it matters only for the continuous (flow) case, where the WBA is
carried as an extra ODE `dW/dζ = g(ζ/T) · h`.

## Cross-validation: two unrelated tools locate the same island

At `K = 6.9` in the Chirikov–Taylor standard map there is a stable accelerator
fixed point. WBA classifies its interior as regular (`median dig = 10.4`) and
measures the island's area on a 500×500 grid. The transport axis (**D1**, in
sv-2/`d1_standard_map`, commit `4bc1fbd`) locates the same island independently,
by measuring the diffusion exponent `μ(K)` and asking why `μ → 2` (ballistic
transport).

| | Position θ\* | Trace | Island area / torus |
|---|---|---|---|
| This work (WBA convergence) | 1.9968 | −0.8516 | **0.405 %** |
| D1 = sv-2/`d1_standard_map` (transport, diffusion exponent `μ`) | 1.997 | −0.85 | **0.37 %** |

The position and trace are analytic properties of the accelerator fixed point
(`sin θ* = 2π/K`, `Tr = 2 + K cos θ*`), so both approaches necessarily share
them. The independent quantity is the **area**: 0.37 % from transport scaling
versus 0.405 % from convergence classification, agreeing to within grid
resolution.

The two methods share no common principle — one integrates the variance of
transport, the other integrates a weighted average and counts digits — and they
agree. This was the first genuine cross-check that the toolkit measures a real
object rather than an artifact of one method.

A related sanity check falls out of the same computation: sampling initial
conditions on the `θ = 0` symmetry line misses the accelerator island entirely
(0 hits in 2000 points), while the 2D grid resolves it. A one-dimensional scan
can therefore report the wrong dynamics; the fraction estimates here use a 2D
grid throughout.

WBA and the Benettin (finite-time Lyapunov) classifier agree on 99.99 % of
14 400 initial conditions at `K = 6.9`. At equal iteration count the two cost
about the same here (`t_WBA = 2.12 s` vs `t_Lyap = 2.63 s`); WBA's advantage is
qualitative — machine-precision digits on regular orbits and a crisp cut, rather
than a slow approach of `λ → 0`.

## What the observable must be

The observable fed to WBA is not innocuous. On the torus the classifier uses
`h = sin θ`, a bounded, smooth function of the orbit. The raw action `p` fails
on accelerator modes, where `p` advances by `2π` per step so the two windows
differ by `~2π T` and a regular island is misread as chaotic. The internal
rotation number is instead recovered as the WBA of the wrapped angle advance.

The general lesson is to prefer position-like observables over angle increments.
Differencing angles introduces two hazards independent of this repository:
Nyquist aliasing once the per-step angle advance exceeds the sampling rate, and
noise amplification, since a difference is a derivative. The 3D application
(sv-2) is where this matters in practice; see that repository for the concrete
case.

## What was not reproduced

Duignan & Meiss report a critical perturbation `ε_cr = 0.665` for the
Section 5.3 magnetic-field model (their Chirikov overlap estimate is 1.0).
**This value was not reproduced here.**

The model was implemented as specified — resonances `(m,n)` at `ψ = n/m` with
amplitudes `(ε/21600)·(72,27,25,96,25,27,72)`, and the WBA carried as an ODE
along the field line. The obstacle is that the paper's data do not fix the
amplitude-to-Hamiltonian normalisation (the factor relating the given weights to
the coefficient of `cos(mθ − nζ)`, i.e. the island-width convention). The
benchmark field-line escape threshold depends strongly on that factor:

| Amplitude convention | ε_cr (benchmark escape) |
|---|---|
| literal (coeff = ε/21600·AW) | ≈ 0.426 |
| overlap-calibrated (all pairs overlap at ε = 1) | ≳ 1.2 |
| paper | 0.665 |

The paper's value is bracketed by the two conventions but not pinned. The
*mechanism* does reproduce: as `ε` increases, the `dig` of the confined
benchmark orbit collapses from ≈ 8–9 (regular, trapped) to ≈ 1 (chaotic,
escaping), which is WBA detecting the confining KAM barrier break. Only the
numeric `ε` is convention-dependent.

**No tuning was performed to force agreement.** The discrepancy is recorded as
it stands, and pinning `0.665` is left to the paper's exact field-line
Hamiltonian. Details in `LIMITATIONS.md` §6 and `future_questions.md`.

## Repository layout

| Path | Contents |
|---|---|
| `wba/weight.py` | Bump weight `g_w`, normalisation constant `C` |
| `wba/birkhoff.py` | Discrete WBA weighted sum, two-window `dig_T` |
| `wba/standard_map.py` | Standard map, WBA sweep (`h = sin θ`), Benettin Lyapunov |
| `wba/magnetic_field.py` | Section 5.3 field-line Hamiltonian, WBA as an ODE |
| `scripts/01_standard_map.py` | `dig` maps, chaos fractions, WBA vs Lyapunov |
| `scripts/02_accelerator_island.py` | K=6.9 island, area, line-vs-grid |
| `scripts/03_magnetic_field.py` | ε_cr scan, chaos-fraction (Fig-14) analog |
| `results/` | JSON summaries and figures |
| `RECORD_B2_WBA.md` | Full results, numbers, pass/fail |
| `LIMITATIONS.md` | Caveats — read first |
| `future_questions.md` | Out-of-scope items, quarantined |

Reproduce with:

```
pip install numpy scipy matplotlib
python3 wba/weight.py                     # C ≈ 142.2503758, ∫g = 1
python3 scripts/01_standard_map.py        # dig maps, chaos fractions, WBA vs Lyapunov
python3 scripts/02_accelerator_island.py  # K=6.9 island, 0.405 %, line-vs-grid
python3 scripts/03_magnetic_field.py      # ε_cr scan, Fig-14 analog
```

## Downstream use

WBA yields the rotation number `ι` as a byproduct of the same integration used
for the classifier. It was applied to 3D QUASR vacuum fields in sv-2 (including
rotation-number profiles), and alongside converse-KAM to Paul's critical-overlap
fields in sv-4.

## Related repositories

- [sv1_comp](https://github.com/kmh21012799-cpu/sv1_comp) — full narrative
- [sv-1](https://github.com/kmh21012799-cpu/sv-1) — standard-map diffusion **coefficient** `D(K)` (a separate transport study; not the `μ` study cross-checked above)
- [sv-2](https://github.com/kmh21012799-cpu/sv-2) — QUASR vacuum configurations (contains `d1_standard_map` = **D1**, the diffusion-**exponent** `μ(K)` study cross-checked above)
- [sv-4](https://github.com/kmh21012799-cpu/sv-4) — converse-KAM 3D + V_PD

## References

- Duignan, N. & Meiss, J. D. (2023) *Physica D* **449**, 133749.
- Sander, E. & Meiss, J. D. (2020) *Physica D* **411**, 132569. *(background on Birkhoff averages for area-preserving maps; listed for context, not cited by the code in this repository.)*
