# LIMITATIONS — B2 WBA (read first)

This is a **reproduction-only** exercise. No novelty, no discovery is claimed.
The goal was to stand up a working Weighted Birkhoff Average (WBA) tool and
reproduce benchmarks from Duignan & Meiss (2023, Physica D 449, 133749).

## Limitations the authors themselves state (know before trusting a result)

1. **"Chaotic ⇒ slow convergence" is a heuristic, not a theorem.** Only the
   *super-convergence of regular (Diophantine) orbits* is proved. A small
   `dig_T` is strong evidence of chaos but not a proof. We use it as a
   classifier, not a certificate.
2. **The observable `h` is not uniquely determined.** D&M report that using an
   approximate invariant as `h` did *not* sharpen the contrast (a negative
   result). We use simple observables (`h = sin θ` for the standard map,
   `h = ψ` for the field) and do not tune them.
3. **Super-convergence needs a Diophantine rotation vector.** Near rationals and
   low-order resonances it degrades; boundary orbits are rare but real.
4. **The authors call their chaos-fraction indicator "crude."** Our chaos
   fractions inherit that; they depend on the `dig_T < 5` cut, on `T`, and on
   the grid. Treat them as coarse, resolution-dependent estimates.

## Limitations specific to this reproduction

5. **No D1 assets were present in this repository.** The task refers to an
   existing D1 standard-map / Lyapunov codebase and to D1's measured
   accelerator-island fraction (0.37%) and μ-scaling. Those files are **not in
   this repo** (it contained only `README.md`). Everything here — standard map,
   Benettin Lyapunov, field-line flow — was implemented from scratch. The
   "0.37%" and the μ=1-vs-μ=2 symmetry-line facts are treated as **external
   claims from the task/D1** and compared against, not re-derived from D1 code.
6. **`ε_cr = 0.665` for the §5.3 magnetic field was NOT uniquely reproduced.**
   The task gives resonance amplitudes as `(ε/21600)·(72,27,25,96,25,27,72)`
   but does **not** fix the amplitude→Hamiltonian normalisation (the factor
   relating the given weights to the coefficient of `cos(mθ−nζ)` in `H`, i.e.
   the island-width convention). The benchmark single-field-line escape
   threshold is strongly sensitive to that factor:
   * literal reading (coeff = ε/21600·AW): escape onset ≈ 0.42
   * overlap-calibrated reading (coeff set so all pairs Chirikov-overlap at
     ε=1): escape onset ≳ 1.2 (island-dominated, non-monotonic)
   The paper's 0.665 is **bracketed** by these conventions but not pinned. The
   *mechanism* (WBA detects the confining-KAM-barrier break via the `dig`
   collapse of the benchmark orbit) is reproduced; the specific number is not,
   pending the paper's exact field-line Hamiltonian. This is logged as an
   **unreproduced anchor**, honestly, rather than tuned to hit 0.665.
7. **`w = 1` is fixed throughout.** No weight tuning. Reported spikes (if any)
   are recorded, not optimised away.
8. The standard map is treated on the 2-torus (θ, p mod 2π). On the torus the
   K>2π accelerator fixed points are ordinary elliptic islands; the "ballistic
   / accelerator" character is a cylinder statement. WBA flags the *island*
   (torus object); the transport consequence (μ=2) is the D1 cylinder result.
