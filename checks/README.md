# Checks (Python)

Independent linear calculations used alongside the Julia scripts. Equations follow
Baker & Mashayek (2021, J. Fluid Mech. 924, A17, doi:10.1017/jfm.2021.627); the code is my own.

| Script | What it computes | Result used in validation_log.md |
|---|---|---|
| `lid_sponge_steady.py` | Steady linear response with a rigid lid and the simulation's sponge (200 m, 600 s), f = 0 | drag/ref 0.915 (sponge as simulated), 1.014 (600 m sponge), about 0 with no damping |
| `impulsive_start_f0.py` | Time-marched linear response to an impulsive start, f = 0, tall domain, mean over t = 4-6 h | drag/ref versus height (the `linear` list in flux_profile.jl) |
| `impulsive_start_lid_sponge.py` | Same, f = 0, but in the simulation's own domain: rigid lid at 1000 m plus the top sponge (200 m, 600 s) | drag/ref at 75 m above the crest, t = 4-6 h: 0.968 |
| `impulsive_start_rotating.py` | Same with rotation f = 1.4e-4, f = 0 for comparison | rot/non-rot ratio 0.69 at 75 m above the crest (the `linear_rot` list) |

Run time: under a few minutes each. Parameters: N = 2e-3, U = 0.2, h0 = 20, sigma = 500.
