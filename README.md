# Lee-wave benchmark in Oceananigans.jl

A small two-dimensional benchmark of internal (lee) wave generation by an isolated ridge in
uniformly stratified flow, with and without rotation. A nonlinear Oceananigans simulation is
compared with linear lee-wave theory.

**Status: work in progress.** The simulated wave structure, the way the flux changes with height
and the effect of rotation agree with linear theory. The simulated wave drag does **not** agree in
magnitude: it is about 1.3 times the linear prediction on every grid I could afford, and I have not
yet explained why (see [Results](#results) and `validation_log.md`).

Author: Ilemona Sunday Omeje.

## What is simulated

A steady current `U` flows over a Gaussian ridge `h(x) = h0 exp(-x^2 / 2 sigma^2)` on the floor of a
2-D (x, z) domain. Stratification is uniform (buoyancy frequency `N`). The ridge radiates internal
waves that carry momentum upward; the ridge feels this as a wave drag.

| Quantity | Value |
|---|---|
| Background flow U | 0.2 m/s |
| Buoyancy frequency N | 2e-3 1/s |
| Ridge height h0, half-width sigma | 20 m, 500 m |
| Lee-wave Froude number N h0 / U | 0.2 |
| Domain | 20 km x 1 km, periodic in x, rigid lid |
| Coriolis parameter f | 0 and 1.4e-4 1/s (about 80 N) |
| Absorbing layers | 2 km at each x-edge and the top 200 m; relaxation time 600 s on u (to U), v, w and the buoyancy perturbation |
| Numerics | Oceananigans `NonhydrostaticModel`, WENO advection, RK3, time step 20 s, CG pressure solver, staircase (grid-fitted) immersed ridge |
| Rotating case | f-plane plus the geostrophic force `f U` in the v equation, so the mean current stays steady |

The diagnostic is the conserved Eliassen-Palm flux `F = <u'w'> - f <v'b'> / N^2`, integrated over
`|x| < 8 km` (outside the absorbing layers) at 50-100 m above the crest and averaged over the last
2 h of a 6 h run. The drag is `-F`. It is normalised by `N U h0^2`.

## Theory

`theory.jl` implements steady linear theory (Bell 1975) in the form given by Baker & Mashayek (2021),
as my own code (not their solver). It uses their equations (2.8), (2.10), (2.18), (2.19), (2.21)-(2.23),
(2.47), (2.49) and (3.1); equation numbers refer to the published paper. For this ridge the linear drag
is 0.980 `N U h0^2` without rotation and 0.748 `N U h0^2` with f = 1.4e-4.

Independent Python checks are in `checks/` (see `checks/README.md`). They include a time-marched linear
model of an impulsive start, in a tall domain and in the simulation's own domain (lid plus sponge), which
is needed because after 4-6 h the waves have not yet reached their steady state, especially with rotation.

## Results

Drag normalised by `N U h0^2`, mean over t = 4-6 h. "Ratio sim/linear" compares with my time-matched
linear prediction at 75 m above the crest. Full details, including runs that failed, are in
`validation_log.md`.

| Case | dx, dz (m) | Drag / ref | Ratio sim/linear |
|---|---|---|---|
| Linear theory, f = 0 (steady) | | 0.980 | |
| f = 0, h0 = 20 m | 39, 3.9 | 1.455 | 1.46 |
| f = 0, h0 = 20 m | 39, 1.95 | 1.366 | 1.37 |
| f = 0, h0 = 20 m | 19.5, 3.9 | 1.313 | 1.32 |
| f = 0, h0 = 20 m | 19.5, 1.95 | 1.317 | 1.32 |
| f = 0, h0 = 40 m (Fr_L = 0.4) | 39, 3.9 | 1.419 | 1.42 |
| Linear theory, f = 1.4e-4 (steady) | | 0.748 | |
| f = 1.4e-4, h0 = 20 m | 19.5, 3.9 | 0.931 | 1.35 |

What the results show:

- **Wave structure and propagation agree with linear theory.** The flux profile with height has the same
  shape as the linear prediction, and the ratio sim/linear is nearly constant from 25 to 300 m above the crest.
- **Rotation behaves as expected.** The drag with rotation divided by the drag without it, on the same
  grid, is 0.709. The time-matched linear prediction is 0.691; the steady value is 0.763 (the simulation
  has not reached steady state at 6 h).
- **The magnitude does not agree.** The simulated drag is about 32 % above the linear prediction. The offset
  falls from 46 % to 32 % as the grid is refined from dx = 39 m to 19.5 m, but halving dz at dx = 19.5 m changes
  nothing, so it is not removed by refinement in the range I could run. It barely depends on ridge
  height (Fr_L 0.2 to 0.4 changes the ratio by about 4 %).
- The 10 % agreement I set as the target at the start is **not met**.

Checked and not the cause (details in `validation_log.md`): the transient after the impulsive start,
reflection from the lid and sponge (about -3 % to -7 %), nonlinearity as the main factor, and the flux
diagnostic itself (applied to the exact linear solution on the staggered grid it returns 0.979 against 0.980,
checked in Python).

Not yet tested: the advection scheme, the time step, and how the immersed ridge is represented. Whether the
offset is a model or numerical effect, or a difference between the simulation and the linear problem that I
have not identified, is open.

## What this does not show

This is a numerical verification exercise. It does not validate the wave energy flux, test
non-uniform background flow, the surface reflection that is the subject of Baker & Mashayek (2021),
or any parameterisation.

## Running it

Tested with Julia 1.13.1, Oceananigans 0.113.5, CairoMakie 0.15.15, JLD2 0.6.7 and FFTW 1.10.0 (pinned in
`Project.toml` and `Manifest.toml`).

```
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project theory.jl                          # linear theory, a few seconds
julia --project run_case.jl 0.0 test               # 30-minute test run, about 1 min
julia --project run_case.jl 0.0 nx=1024            # f = 0, 6 h, about 16 min
julia --project run_case.jl 1.4e-4 nx=1024         # with rotation, about 15 min
julia --project compare.jl lee_f0_nx1024.jld2 0.0
julia --project flux_profile.jl lee_f0_nx1024.jld2 0.0
```

At the Julia prompt, wrap each command in `run(`...`)`. The second argument of `compare.jl` and
`flux_profile.jl` must be the same f used for the run. Run times are for a laptop with an Intel i5-8250U
and 8 GB of memory. Options of `run_case.jl`: `test`, `centered`, `nx=...`, `nz=...`, `h0=...`. Output
files (`lee_*.jld2`, up to about 100 MB each) are not kept in the repository.

| File | Purpose |
|---|---|
| `params.jl` | All parameters in one place; `describe(p)` prints derived checks |
| `theory.jl` | Linear theory (the "answer key") |
| `run_case.jl` | The Oceananigans simulation |
| `compare.jl` | Simulation against theory: table and two figures |
| `flux_profile.jl` | Flux against height, next to the linear prediction |
| `validation_log.md` | Every run, including failures, and the open question |
| `CHANGELOG.md` | How the scripts evolved from the starting script |
| `figures/` | Comparison figures |
| `checks/` | Independent Python calculations |

## References

- Bell, T. H. (1975). Topographically generated internal waves in the open ocean. *J. Geophys. Res.* 80, 320-327.
- Baker, L. E. & Mashayek, A. (2021). Surface reflection of bottom generated oceanic lee waves.
  *J. Fluid Mech.* 924, A17. doi:10.1017/jfm.2021.627. Their solver is at
  https://github.com/loisbaker/lee-wave-solver (not used here).
- Mayer, F. T. & Fringer, O. B. (2017). An unambiguous definition of the Froude number for lee waves in the
  deep ocean. *J. Fluid Mech.* 831, R3.
- Ramadhan, A. et al. (2020). Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs.
  *J. Open Source Softw.* 5(53).

## Acknowledgement of assistance

I developed the code, theory checks and documentation with the help of an AI assistant (Claude,
Anthropic). I ran all the simulations and am responsible for the content.
