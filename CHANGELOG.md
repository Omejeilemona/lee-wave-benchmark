# Changelog

Changes from the starting script (`archive/internal_wave_demo.jl`, kept unmodified) to the current
scripts. It was added to the repository after the later versions, so the git history itself begins at the
theory stage. `internal_wave_demo_corrected.jl` is an updated copy of the starting script with only its
parameters corrected (changes marked in the file; the full 18 h run gave about twice the linear drag, see `validation_log.md`). The first commits hold the first versions of the
theory and simulation scripts; edits made before them are listed below, **reconstructed from
development notes, not from commits.** See `archive/about_starting_script.md` for what was wrong
with the starting script.

## Physics / set-up
- Stratification N: 1e-2 -> 2e-3 s^-1, so the lee-wave Froude number Fr_L = N*h0/U goes
  from 1.0 to 0.2 (linear theory needs Fr_L << 1; Baker & Mashayek 2021).
- Time step 5 s -> 20 s and run length 18 h -> 6 h (U*dt/dx = 0.10).
- Sponge now also relaxes the buoyancy perturbation b to 0 (stops buoyancy anomalies
  re-entering through the periodic edge).
- Rotating case: added the geostrophic-balance force f*U0 to the v equation. Without it
  the Coriolis force rotates the whole current (inertial oscillation, period 2*pi/f).
- Advection: Centered() -> WENO(); `centered` flag keeps the old option.

## Diagnostics
- Reference solution: `theory.jl` (linear theory, with and without rotation) replaces the
  single formula N*U*h0^2. For rotation the conserved flux F = u'w' - f v'b'/N^2 is used.
- Output: u, v, w, b every 15 min (v and b are needed for F); GIF removed.
- Flux measured 50-100 m above the ridge crest, over the interior (outside the sponges),
  averaged over the last 2 h, instead of one level at z = -600 m.
- `compare.jl`: simulation vs theory table and two figures.
- `flux_profile.jl`: flux versus height, next to a linear impulsive-start prediction.

## Code structure
- `params.jl` holds every number; `describe(p)` prints derived checks.
- `run_case.jl` options: `test`, `centered`, `nx=...`, `nz=...`, `h0=...`. File names carry the
  resolution and ridge height so runs do not overwrite each other.
- `compare.jl` and `flux_profile.jl` read the grid size from the output file and the ridge
  height from the file name.
- `halo = (4, 4)` set on the grid (y is Flat, so two entries).
