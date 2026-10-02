# The starting script

`internal_wave_demo.jl` is the script I started from, kept here exactly as it was. It is **not** used for any
result in this repository, and I have not run it in the pinned environment. It was added to the repository
after the later scripts, so the git history starts at the theory stage (see `CHANGELOG.md`).

## Why it was replaced

- **Not in the linear regime.** With N = 1e-2 1/s, U = 0.2 m/s and h0 = 20 m, the lee-wave Froude
  number N h0 / U is 1.0, not the 0.3 stated in its comment. Linear theory needs a value well below 1.
- **Wrong reference for the rotating case.** Rotation (f = 1.4e-4) was on, but the comparison used the
  non-rotating hydrostatic formula N U h0^2. With rotation, u'w' alone is not conserved with height; the
  conserved quantity is F = u'w' - f v'b'/N^2 (Baker & Mashayek 2021, section 2.8).
- **No balance for the Coriolis force on the mean current.** With f not zero the whole current would rotate
  (an inertial oscillation) unless a balancing force is added.
- **One measurement level.** The flux was taken at a single level (z = -600 m), integrated over the whole
  width including the absorbing layers, after 18 h.
- **Absorbing layers did not act on buoyancy**, so buoyancy anomalies could re-enter through the periodic edge.
- **Numerics.** Centered advection with no closure; a time step and run length longer than needed.
- Minor: it clears the terminal at start, and it carries some duplicate forcing methods.

The changes are listed in `CHANGELOG.md`.

An updated copy with only the parameters corrected is `../internal_wave_demo_corrected.jl`; its header lists the
changes and what was left as it was.
