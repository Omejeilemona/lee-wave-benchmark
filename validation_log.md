# Validation log

Reference: N*U*h0^2 = 0.16 m^3/s^2 (N = 2e-3, U = 0.2, h0 = 20, sigma = 500).
Linear theory (theory.jl): drag/ref = 0.980 for f = 0; 0.748 for f = 1.4e-4 using the
conserved flux F = u'w' - f v'b'/N^2 (u'w' alone would give 1.070).
Tolerance chosen for this demo: about 10 %. Dates: fill in.

| # | Date | What | Command | Result: drag/ref | Notes |
|---|------|------|---------|------------------|-------|
| 0 | | Starting script `internal_wave_demo.jl` | -- | not compared | N = 1e-2 gives Fr_L = 1.0; compared a rotating run with a non-rotating formula |
| 1 | | Theory self-test (Julia) | `julia --project theory.jl` | 0.980 (f=0); 0.748 (f=1.4e-4) | w(x,0) = U dh/dx to 1e-14. Matches an independent Python port. |
| 2 | | f=0, 30-min timing test | `run_case.jl 0.0 test` | 0.064 | Only 3 snapshots; waves had not reached the measurement band. Not a result. |
| 3 | | f=0, 6 h, Nz=256, WENO | `run_case.jl 0.0` | **1.455 (+48.5 %)** | Plateau from 4 h, drift -1.5 %. 7.6 min. |
| 4 | | f=0, 6 h, Nz=512, WENO | `run_case.jl 0.0 nz=512` | **1.366 (+39.4 %)** | Drift -1.6 %. 17.8 min. |
| 5 | | Flux vs height, runs 3 and 4 | `flux_profile.jl` | ratio sim/linear = 1.46 (Nz=256), 1.37 (Nz=512), constant from 25 to 500 m | Wave propagation matches; amplitude at generation is about 17-21 % high. |
| 6 | | f=0, Nz=512, h0=40 m | `run_case.jl 0.0 nz=512 h0=40` | aborted | Stalled after about 50 steps; max\|w\| = 8.7e-2 m/s (vertical Courant number about 0.9, still growing). Not a result. |
| 7 | | f=0, 6 h, Nz=256, h0=40 m (ridge 10 cells, Fr_L = 0.4) | `run_case.jl 0.0 nz=256 h0=40` | **1.419 (+44.9 %)** | Drift +4.7 %. Profile ratio sim/linear 1.42-1.43 (75-200 m). 8.8 min. Same cells per ridge height as run 4 (1.366): doubling Fr_L raised the ratio only about 4 %, so nonlinearity is a small part of the excess. |
| 8 | | f=0, 6 h, Nx=1024, Nz=256, h0=20 m | `run_case.jl 0.0 nx=1024` | **1.313 (+34.1 %)** | Drift -1.6 %. Profile ratio sim/linear 1.32-1.34 (25-300 m), 1.28 at 400-500 m. 16 min. Halving dx lowered the excess more (ratio sim/linear 1.46 -> 1.32) than halving dz (1.46 -> 1.37). Fitting excess = a + b*dx + c*dz through runs 3, 4, 8 gives a = 0.00; with dx^2, dz^2 instead a = 0.15. Both predict about 1.23 at dx = 19.5 m, dz = 1.95 m. Three points fit three parameters exactly, so this is a plausibility argument, not proof. |
| 9 | | f=1.4e-4, 6 h, Nx=1024, Nz=256, h0=20 m | `run_case.jl 1.4e-4 nx=1024` | **0.931 (+24.5 % vs steady theory 0.748)** | Drift +2.9 %. 14.7 min. Drag(rot)/drag(f=0) on the same grid: 0.931/1.313 = 0.709; linear impulsive-start prediction 0.691 (steady theory 0.763). Profile ratio sim/linear 1.40 (25 m) to 1.32 (300 m); 1.35 at 75 m against 1.32 for f=0. u'w' alone gives 1.254 against F = 0.931; the f v'b'/N^2 term is 0.32 in both simulation and steady theory (units of N U h0^2). Check fig_flux_lee_frot_nx1024.png for any oscillation with the 12.5 h inertial period. |
| 10 | | f=0, 6 h, Nx=1024, Nz=512, h0=20 m | `run_case.jl 0.0 nx=1024 nz=512` | pending | Balanced refinement. Predicted ratio sim/linear about 1.23 if the errors are additive discretisation errors. About 40 min. |

## Checks (own Python calculations, in `checks/`)
- Steady linear model with rigid lid and the same sponge (200 m, 600 s): drag/ref = 0.915;
  600 m sponge: 1.014; lid with no damping: about 0 (as in Baker & Mashayek 2021, sec 2.8).
  So lid/sponge reflection explains at most about -7 %, not +39 %.
- Linear impulsive-start model (no sponge, tall domain): drag at 75 m above the crest reaches
  about 0.99 by 5 h with no overshoot.
- Same model with rotation (f = 1.4e-4), mean over t = 4-6 h, 75 m above the crest: drag/ref = 0.689
  against 0.997 for f = 0, a ratio of 0.69 (steady theory: 0.763). Low-wavenumber waves have not
  equilibrated by 4-6 h.

## Open question
Why the simulated flux exceeds linear theory by 30-45 % on the grids used. Runs 3-8: the excess barely
depends on ridge height (Fr_L 0.2 -> 0.4: +4 %) but falls with both horizontal and vertical resolution
(ratio sim/linear 1.46 -> 1.37 for dz halved, 1.46 -> 1.32 for dx halved). Lid/sponge reflection explains
at most about -7 %. Current reading: mostly discretisation error of the staircase ridge, converging towards
linear theory; the converged residual is uncertain (0 to about 15 %). Run 10 tests this.
