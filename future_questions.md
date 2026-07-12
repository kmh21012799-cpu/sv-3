# future_questions.md — quarantined, out of B2 scope

Items that came up while implementing WBA but are explicitly **out of scope**
for B2. Recorded here so they don't leak into the reproduction. Not to be acted
on without a separate decision.

- **converse-KAM overlay.** The sharpest use of WBA is where it disagrees with
  converse-KAM: WBA marks an island "regular/alive", converse-KAM marks it
  "no barrier/removed". The split locates the island. B2 forbids this overlap
  ("이게 진짜 균열이지만 지금 아니다"). Revisit only after B2 passes.

- **`ε_cr` amplitude-normalisation ambiguity.** The §5.3 field amplitudes are
  given as weights `/21600` but the map from weights to the `cos` coefficient
  (island-width convention) is not fixed by the task; `ε_cr` moved from ~0.42
  (literal) to ~1.2 (overlap-calibrated). Resolving this needs the exact D&M
  field-line Hamiltonian. Worth pinning down — but as a documentation fetch,
  not a tuning exercise.

- **`ε_cr` as a proper barrier (noble KAM) destruction, not single-IC escape.**
  The single-field-line escape is island-contaminated and non-monotonic in ε.
  A Greene-residue or WBA-based converse-KAM barrier measure would give a
  cleaner, monotonic `ε_cr`. This edges into converse-KAM territory (above) —
  quarantined.

- **V_PD transport axis.** The §5.3 field is the same as the Paul–Hudson–
  Helander V_PD field. Tempting to bolt on transport. Forbidden in B2.

- **QUASR fields / WBA→μ prediction.** B4 work; the ultimate goal (predict the
  transport exponent μ from island area found by WBA) is explicitly not B2.

- **Weight spikes vs `w`.** Certain `w` produce local spikes when the weighted
  window width matches an integer multiple of the rotation number. Only to be
  *recorded*, never used to pick an "optimal" w. If seen, log location and move
  on.

- **Observable choice.** `h = sin θ` (map) and `h = ψ` (field) were adequate;
  D&M's negative result on approximate-invariant observables suggests limited
  upside, but a systematic `h` sweep was not done.
