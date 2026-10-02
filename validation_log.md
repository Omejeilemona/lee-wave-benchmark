# Validation log

Reference: N*U*h0^2 = 0.16 m^3/s^2 (N = 2e-3, U = 0.2, h0 = 20, sigma = 500).
Linear theory (theory.jl): drag/ref = 0.980 for f = 0; 0.748 for f = 1.4e-4 using the
conserved flux F = u'w' - f v'b'/N^2 (u'w' alone would give 1.070).
Tolerance chosen for this demo: about 10 %. Dates: fill in.

**Status:** diagnosis of the offset is paused after run 10. Rows 11-13 are planned, not run.

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
| 10 | | f=0, 6 h, Nx=1024, Nz=512, h0=20 m | `run_case.jl 0.0 nx=1024 nz=512` | **1.317 (+34.4 %)** | Drift -1.4 %. Profile ratio sim/linear 1.32 (25-150 m), 1.28-1.34 elsewhere. 40.5 min. **Prediction of about 1.23 failed:** halving dz at dx = 19.5 m changed nothing (run 8: 1.313). The convergence extrapolation in row 8 is withdrawn. |
| 11 | | f=0, Nx=1024, Nz=256, Centered advection | `run_case.jl 0.0 nx=1024 centered` | not run | Planned: does the advection scheme matter? Same grid as run 8. |
| 12 | | f=0, Nx=1024, Nz=256, dt = 10 s | (option not yet in the script) | not run | Planned: does the time step matter? |
| 13 | | f=0, Nx=1024, Nz=256, partial-cell ridge | (option not yet in the script) | not run | Planned: does the representation of the ridge matter? Partial cells may not be supported by the pressure solver. |

## Checks (own Python calculations, in `checks/`)
- Steady linear model with rigid lid and the same sponge (200 m, 600 s): drag/ref = 0.915;
  600 m sponge: 1.014; lid with no damping: about 0 (as in Baker & Mashayek 2021, sec 2.8).
  So lid/sponge reflection explains at most about -7 %, not +39 %.
- Linear impulsive-start model (no sponge, tall domain): drag at 75 m above the crest reaches
  about 0.99 by 5 h with no overshoot.
- Same model in the simulation's own domain (rigid lid at 1000 m + the top sponge, f = 0), mean over
  t = 4-6 h: drag/ref at 75 m above the crest = 0.968 (0.997 in the tall domain). So the lid, the sponge and
  the 1000 m depth do not explain the +30 % in the transient either.
- Same model with rotation (f = 1.4e-4), mean over t = 4-6 h, 75 m above the crest: drag/ref = 0.689
  against 0.997 for f = 0, a ratio of 0.69 (steady theory: 0.763). Low-wavenumber waves have not
  equilibrated by 4-6 h.

## Open question
Why the simulated flux is about 1.32 times linear theory (ratio sim/linear, 75 m above the crest) on every grid
with dx <= 19.5 m and dz <= 3.9 m. Ratio by grid (dx, dz): (39, 3.9) 1.46; (39, 1.95) 1.37; (19.5, 3.9) 1.32;
(19.5, 1.95) 1.32. So the excess is not removed by refinement in this range. It also barely changes with ridge height
(Fr_L 0.2 -> 0.4: +4 %), and the profile with height has the same shape as the linear prediction.
Ruled out so far: nonlinearity as the main cause, spin-up transient, lid/sponge/depth (steady -7 %, transient -3 %),
the flux diagnostic (the same operations applied to the exact linear solution on the staggered grid give 0.979
against 0.980; checked in Python). Rotation behaves as expected (rot/non-rot 0.709 against 0.691).
Remaining candidates (planned tests: rows 11-13): advection scheme, time step, treatment of the immersed ridge, a model setting,
or a conceptual difference between the simulation and the linear problem that I have not identified.
