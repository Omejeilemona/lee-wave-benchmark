# =====================================================================
# params.jl  --  every number used by the benchmark lives in this file
# =====================================================================
# Used by: theory.jl (the "answer key"), run_case.jl (the simulation),
#          compare.jl (simulation vs theory).
# Units: SI throughout (metres, seconds).
#
# The set-up: a steady current U flows over a small Gaussian ridge in
# water whose density increases with depth (buoyancy frequency N).
# Gravity pulls lifted water back down, so the ridge radiates a steady
# pattern of internal "lee" waves that carry momentum upward. The ridge
# feels this as a drag. Rotation (Coriolis parameter f) is run as a
# second case.
#
# Where the numbers come from (be ready to say this in an interview):
#   U, h0, sigma : taken from the original script.
#   f            : 0 for the validation case. 1.4e-4 s^-1 is roughly
#                  2*Omega*sin(80 deg) = 1.4365e-4 (Omega = 7.2921e-5),
#                  the value used in the original script.
#   N            : 2e-3 s^-1 chosen so that the lee-wave Froude number
#                  Fr_L = N*h0/U = 0.2 (linear regime). The original
#                  N = 1e-2 gave Fr_L = 1.0, which is not small.
#                  (Fr_L << 1 for linear theory: Baker & Mashayek 2021.)
#   grid, dt,... : chosen to keep the run light on an 8 GB laptop.
# =====================================================================

using Printf

const f_arctic = 1.4e-4   # s^-1

"""
    benchmark_params(; f = 0.0, alpha = 1)

Return all parameters as a NamedTuple.
  f     : Coriolis parameter (0.0 = no rotation, f_arctic = rotating case)
  alpha : 1 = nonhydrostatic (matches the simulation), 0 = hydrostatic
"""
function benchmark_params(; f = 0.0, α = 1)
    return (;
        # --- physics ---
        f, α,
        U  = 0.2,        # m/s   background flow
        N  = 2e-3,       # 1/s   buoyancy frequency (uniform stratification)
        h0 = 20.0,       # m     ridge height
        σ  = 500.0,      # m     ridge half-width (Gaussian)
        # --- domain / grid (2D: x and z) ---
        Lx = 20e3,       # m
        Lz = 1000.0,     # m
        Nx = 512,
        Nz = 256,
        # --- time stepping and output ---
        Δt = 20.0,                 # s
        stop_time = 6 * 3600.0,    # s  (6 hours)
        output_interval = 900.0,   # s  (15 minutes)
        # --- sponge layer (absorbs waves near the top and the x-edges) ---
        sponge_width     = 2e3,    # m   at each x-edge
        top_sponge       = 200.0,  # m   thickness at the top
        sponge_timescale = 600.0,  # s   relaxation time
        # --- where/when to measure the flux in compare.jl ---
        avg_window    = 2 * 3600.0,    # s   average over the last 2 hours
        z_above_crest = (50.0, 100.0), # m   measure 50-100 m above the crest
    )
end

"""
    describe(p)

Print derived quantities so you can see at a glance whether the set-up
makes sense (linear regime? ridge resolved? run long enough?).
"""
function describe(p)
    dx = p.Lx / p.Nx
    dz = p.Lz / p.Nz
    Fr = p.N * p.h0 / p.U
    k  = 1 / p.σ                              # typical ridge wavenumber
    cgz = k * p.U^2 / p.N                     # vertical group velocity, hydrostatic non-rotating
    println("---------------- benchmark parameters ----------------")
    @printf("f = %.2e 1/s   (%s)\n", p.f, p.f == 0 ? "no rotation" : "rotating")
    @printf("Lee-wave Froude number  Fr_L = N*h0/U = %.2f   (want << 1)\n", Fr)
    @printf("Grid: dx = %.1f m, dz = %.2f m;  ridge height = %.1f grid cells\n", dx, dz, p.h0 / dz)
    @printf("Hydrostatic check  N*sigma/U = %.1f   (want >> 1)\n", p.N * p.σ / p.U)
    @printf("Rotation cutoff  f/U = %.2e 1/m  ->  (f/U)*sigma = %.2f\n", p.f / p.U, p.f / p.U * p.σ)
    @printf("Vertical wavelength 2*pi*U/N = %.0f m   (domain depth %.0f m)\n", 2π * p.U / p.N, p.Lz)
    @printf("Vertical group velocity (k = 1/sigma, f = 0) = %.1e m/s -> %.2f h to climb 75 m\n",
            cgz, 75 / cgz / 3600)
    @printf("Advective CFL  U*dt/dx = %.2f ;  steps per run = %d\n",
            p.U * p.Δt / dx, round(Int, p.stop_time / p.Δt))
    Fr > 0.3 && @warn "Fr_L > 0.3: linear theory may not apply."
    p.h0 / dz < 4 && @warn "Ridge is under 4 grid cells tall: poorly resolved."
    println("------------------------------------------------------")
end
