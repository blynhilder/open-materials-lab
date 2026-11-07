# Open Materials Lab — Progress

> Teaching simulation / reproducible crystallography computations.
> No fabricated data. No novelty claims. PRs for review; never push directly to `main`.

## Increment 1 (done) — Crystal structure → powder XRD simulation

**Branch / PR:** `increment-01-xrd` → `main` (open PR #1); local commit `3eeb860` on `master`.

- From-scratch powder XRD (`src/xtal/`, numpy-first) for Si (diamond), Al (FCC), α-Fe (BCC).
- COD CIFs: 9011998 / 9012002 / 9008536 with SHA-256 in `data/SOURCES.md`.
- Cross-validation vs pymatgen `XRDCalculator`: max |Δ2θ| ~ 10⁻¹⁴ °; intensity differences < 1/100 from form-factor tables; identical when using the same f(s).
- 53 pytest tests green; figures and validation CSVs under `figures/` and `results/`.
- Limitations: B = 0, no absorption / preferred orientation / Kα2 / instrument function; peak positions not yet compared to experimental standards.

## Increment 2 (this week) — NIST SRM 640 Si lattice / peak-position benchmark

**Branch:** `increment-02-nist-srm640`.

Goal (roadmap step **a**): benchmark against experiment/standards — NIST SRM 640 Si certified lattice parameter, with cited sources.

Delivered:
1. Downloaded NIST SRM 640e & 640f certificates (2026-11-07 Asia/Shanghai); SHA-256 recorded in `data/SOURCES.md` and `data/nist/SRM_640{e,f}.json`.
2. Machine-readable transcriptions of certified *a* and Table 1 peak lists (`data/nist/*.json`, `*_table1.csv`).
3. Helper `cubic_bragg_two_theta` + script `scripts/run_nist_srm640_benchmark.py`.
4. Pytest suite `tests/test_nist_srm640.py` (certified *a*, Table 1 recomputation ≤ 0.001°, COD *a*, PDF checksums).
5. Results: `results/nist_srm640_benchmark.md`, CSVs, figure `figures/Si_nist_srm640e_peak_positions.png`.

### Key validation numbers (from certificates + Bragg recomputation)

| Quantity | Value |
|---|---|
| COD Si a (9011998, 298.15 K) | 5.430941 Å |
| NIST SRM 640e certified a (22.5 °C) | 5.431179 Å ± 0.00008 Å |
| NIST SRM 640f certified a (22.5 °C) | 5.431144 Å ± 0.00008 Å |
| Δa (COD − 640e) | −0.000238 Å (−43.82 ppm) |
| Δa (COD − 640f) | −0.000203 Å (−37.38 ppm) |
| Certificate Table 1 λ | 1.5405929 Å (0.15405929 nm) |
| max \|Δ2θ\| recompute 640e Table 1 from certified a+λ | ~0.00048° (< 0.001°) |
| max \|Δ2θ\| recompute 640f Table 1 from certified a+λ | ~0.00050° (< 0.001°) |

Clear wording: this is a **teaching simulation + literature/COD vs NIST certified standard comparison**. We did **not** measure SRM powder. Table 1 entries are NIST *information values*.

### Limitations (Increment 2)

- No thermal expansion correction between 22.5 °C (NIST) and 298.15 K (COD).
- Only peak *positions* (Bragg geometry) are checked; intensities / Debye–Waller / absorption / instrument function remain out of scope.
- Table 1 information values have no stated uncertainty on the certificate.

## Next step (Increment 3 proposed) — (b) isotropic Debye–Waller B factors

Add isotropic Debye–Waller factors `exp(−B sin²θ / λ²)` with **literature-cited** B values for Si / Al / Fe (or elemental defaults with sources), recompute relative intensities, and document the effect vs B = 0. Keep all sources cited with URL/ID/date; no fabricated B values.

## Roadmap reminder (from Increment 1 plan)

- (a) ~~benchmark vs experiment/standards (NIST SRM 640)~~ ← Increment 2
- (b) isotropic Debye–Waller B factors with literature-cited B ← **next**
- (later) absorption, Kα2, instrument function, non-cubic systems, etc.
