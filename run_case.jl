# =====================================================================
# run_case.jl  --  Oceananigans simulation: steady flow over a Gaussian ridge
# =====================================================================
# Usage from the julia> prompt, inside this folder (one case per run):
#
#   run(`julia --project run_case.jl 0.0 test`)       # f = 0, 30-min TIMING test
#   run(`julia --project run_case.jl 0.0`)            # f = 0, full run (validation case)
#   run(`julia --project run_case.jl 1.4e-4`)         # rotating case
#   run(`julia --project run_case.jl 0.0 centered`)   # optional: Centered advection instead of WENO
#   run(`julia --project run_case.jl 0.0 nz=512`)     # optional: finer vertical grid (file: lee_f0_nz512.jld2)
#   run(`julia --project run_case.jl 0.0 nz=512 h0=40`) # optional: ridge height 40 m (file: lee_f0_nz512_h40.jld2)
#
# Output file names:  lee_f0.jld2 / lee_frot.jld2  (+ _centered, _test if used)
#
# What is simulated (2D: x along the flow, z up):
#   - uniform background flow U over a Gaussian ridge on the bottom
#   - uniform stratification N (background buoyancy N^2 z; the model's
#     tracer b is the PERTURBATION from this background)
#   - optional rotation f, with the matching geostrophic force (see below)
#   - a sponge near the top and at the two x-edges absorbs the waves
#
# Changes from internal_wave_demo.jl are listed in the README / chat log.
# =====================================================================

using Oceananigans
using Oceananigans.Solvers: ConjugateGradientPoissonSolver
using Printf
include(joinpath(@__DIR__, "params.jl"))

# ---- command-line options ---------------------------------------------
f_val     = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.0
flags     = ARGS[2:end]
test_mode = "test" in flags
centered  = "centered" in flags

nz_flag = findfirst(startswith("nz="), flags)
Nz_val  = nz_flag === nothing ? 256 : parse(Int, split(flags[nz_flag], "=")[2])

h0_flag = findfirst(startswith("h0="), flags)
h0_val  = h0_flag === nothing ? 20.0 : parse(Float64, split(flags[h0_flag], "=")[2])

p = benchmark_params(; f = f_val, Nz = Nz_val, h0 = h0_val)
describe(p)

res_tag  = (Nz_val == 256 ? "" : "_nz$(Nz_val)") * (h0_val == 20.0 ? "" : "_h$(round(Int, h0_val))")   # keeps runs in different files
tag      = (p.f == 0 ? "f0" : "frot") * res_tag * (centered ? "_centered" : "") * (test_mode ? "_test" : "")
filename = "lee_" * tag * ".jld2"
stop     = test_mode ? 1800.0 : p.stop_time
@info "Running" f = p.f advection = (centered ? "Centered" : "WENO") stop_time_s = stop output = filename

# ---- 1. grid with immersed Gaussian ridge --------------------------------
underlying_grid = RectilinearGrid(
    size = (p.Nx, p.Nz),
    x = (-p.Lx / 2, p.Lx / 2),
    z = (-p.Lz, 0),
    topology = (Periodic, Flat, Bounded),
    halo = (4, 4),      # (x, z) only: y is Flat. An ImmersedBoundaryGrid needs the extra halo point.
)

h0, σ, Lz = p.h0, p.σ, p.Lz
bottom_topography(x) = -Lz + h0 * exp(-x^2 / (2σ^2))
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_topography))

# ---- 2. stratification (background field) --------------------------------
@inline background_b(x, z, t, params) = params.N2 * z
b_background = BackgroundField(background_b, parameters = (N2 = p.N^2,))

# ---- 3. sponge + geostrophic balance ------------------------------------
# Constants as `const` so the forcing functions are fast and type-stable.
const LX  = p.Lx
const SPW = p.sponge_width
const TOP = p.top_sponge
const U0  = p.U
const TAU = p.sponge_timescale
const GEO = p.f * p.U      # geostrophic balance force, see below

# 0 in the interior, ramping up quadratically to 1 inside the sponge
@inline sponge_mask(x, z) = max(max(0.0, (abs(x) - (LX / 2 - SPW)) / SPW)^2,
                                max(0.0, (z + TOP) / TOP)^2)

# Relax u -> U0, and v, w, b' -> 0 inside the sponge.
@inline u_damping(x, z, t, u) = -sponge_mask(x, z) * (u - U0) / TAU
@inline v_damping(x, z, t, v) =  GEO - sponge_mask(x, z) * v / TAU
@inline w_damping(x, z, t, w) = -sponge_mask(x, z) * w / TAU
@inline b_damping(x, z, t, b) = -sponge_mask(x, z) * b / TAU

# Safety: in case your Oceananigans version passes y to forcing functions on a
# grid with a Flat y dimension, these forward to the versions above.
@inline u_damping(x, y, z, t, u) = u_damping(x, z, t, u)
@inline v_damping(x, y, z, t, v) = v_damping(x, z, t, v)
@inline w_damping(x, y, z, t, w) = w_damping(x, z, t, w)
@inline b_damping(x, y, z, t, b) = b_damping(x, z, t, b)

# GEO = f*U0 is the background pressure-gradient force that balances the
# Coriolis force on the mean flow (geostrophic base state, Baker & Mashayek
# 2021, eq 2.4). WITHOUT it the Coriolis force would rotate the whole
# current (an inertial oscillation, period 2*pi/f = 12.5 h for f = 1.4e-4).
# It is zero when f = 0.
forcing = (
    u = Forcing(u_damping, field_dependencies = :u),
    v = Forcing(v_damping, field_dependencies = :v),
    w = Forcing(w_damping, field_dependencies = :w),
    b = Forcing(b_damping, field_dependencies = :b),
)

# ---- 4. model -----------------------------------------------------------
model = NonhydrostaticModel(
    grid;
    pressure_solver = ConjugateGradientPoissonSolver(grid),
    advection = centered ? Centered() : WENO(),
    timestepper = :RungeKutta3,
    coriolis = p.f == 0 ? nothing : FPlane(f = p.f),
    buoyancy = BuoyancyTracer(),
    tracers = :b,
    background_fields = (b = b_background,),
    forcing = forcing,
)

set!(model, u = p.U)

# ---- 5. simulation and output ------------------------------------------
simulation = Simulation(model, Δt = p.Δt, stop_time = stop)

wall_start = time_ns()
function progress(sim)
    w = sim.model.velocities.w
    @printf("iter %5d   t = %5.2f h   max|w| = %.2e m/s   wall = %.0f s\n",
            iteration(sim), time(sim) / 3600, maximum(abs, interior(w)),
            (time_ns() - wall_start) / 1e9)
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

outputs = merge(model.velocities, model.tracers)     # u, v, w and b
simulation.output_writers[:fields] = JLD2Writer(
    model, outputs;
    filename = filename,
    schedule = TimeInterval(p.output_interval),
    overwrite_files = true,
)

run!(simulation)

wall = (time_ns() - wall_start) / 1e9
@printf("Finished. Wall time = %.0f s (%.1f min) for %.2f simulated hours.\n", wall, wall / 60, stop / 3600)
if test_mode
    @printf("Estimated full 6 h run: about %.0f min.\n", wall / stop * p.stop_time / 60)
end
