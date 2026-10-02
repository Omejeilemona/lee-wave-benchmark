# =====================================================================
# internal_wave_demo_corrected.jl
# The starting script (archive/internal_wave_demo.jl) with its parameters corrected, as a
# single-file demo. Every changed line is marked "CHANGED". The full 18 h run completed (about 1.5 h on
# the author's laptop); see the result at the end of this header.
#
# CHANGED:
#   1. N2 = 1e-4 -> 4e-6 (N = 2e-3 1/s). The original gave N*h0/U = 1.0, not the 0.3 in its
#      comment; now N*h0/U = 0.2, which is in the linear regime.
#   2. f0 = 1.4e-4 -> 0 (no rotation). The reference line (-N U0 h0^2) is the non-rotating formula.
#      Rotation needs a balancing force on the mean current and the conserved flux
#      u'w' - f v'b'/N^2; see run_case.jl and compare.jl for that.
#   3. Time step 5 s -> 20 s (advective CFL 0.10), as in the main benchmark. Run length unchanged.
#   4. Output file names carry "_corrected" so they do not overwrite the original's outputs.
#   5. Legend moved from the bottom right to the top right so that it does not cover the curve.
#
# NOT changed (known limits): Centered advection with no closure; the sponge does not act on buoyancy;
# the flux is taken at one level (z = -600 m) over the whole width; the reference is the hydrostatic
# value (the nonhydrostatic linear value for these parameters is 0.98 of it, from theory.jl).
# RESULT of the full run (validation_log.md, row 14): the simulated flux at z = -600 m crosses the red line
# at about 4 h and averages 2.22 times the red line over the last 2 h (2.3 times the steady linear value of
# 0.98, drift 10 %), with oscillations after about 10 h. I had
# expected 1.4-1.5 times (runs 3-4), so the expectation was wrong: with Centered advection the offset is
# larger than with WENO on the same grid. The flux level, averaging window and duration differ from
# compare.jl, so this is suggestive, not a like-for-like comparison.
# =====================================================================
using Oceananigans
using Oceananigans.Solvers: ConjugateGradientPoissonSolver
using Oceananigans.Units
using CairoMakie
using Printf

# Clear terminal screen
run(Sys.iswindows() ? `cmd /c cls` : `clear`)

# 1. Grid Definition (2D domain: x, z)
Nx, Nz = 512, 256
Lx, Lz = 20kilometers, 1000meters

underlying_grid = RectilinearGrid(
    size = (Nx, Nz),
    x = (-Lx/2, Lx/2),
    z = (-Lz, 0),
    topology = (Periodic, Flat, Bounded)
)

# 2. Topography (Linear Gaussian Ridge)
h0 = 20meters          # Linear regime: N*h0/U = 0.2 << 1 (CHANGED comment: was "~ 0.3")
sigma = 500meters      # Narrow ridge: Uk ~ 4e-4 s⁻¹ > f0

bottom_topography(x) = -Lz + h0 * exp(-x^2 / (2 * sigma^2))

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_topography))

# 3. Arctic Stratification, Coriolis & Flow Parameters
N2 = 4e-6              # CHANGED: N = 2e-3 s⁻¹ (was N2 = 1e-4, N ~ 0.01 s⁻¹, which gives N*h0/U = 1.0)
N = sqrt(N2)
U0 = 0.2               # Background inflow velocity (m/s)
f0 = 0.0               # CHANGED: no rotation (was 1.4e-4, Arctic ~80°N); see header note 2

@printf("Lee-wave Froude number N*h0/U0 = %.2f\n", N * h0 / U0)   # CHANGED: added check

@inline background_b(x, z, t, p) = p.N2 * z
b_field = BackgroundField(background_b, parameters=(N2=N2,))

# 4. Sponge Layer Forcing Functions
top_sponge = 200meters
sponge_width = 2kilometers

@inline function sponge_mask(x, z)
    x_mask = max(0, (abs(x) - (Lx/2 - sponge_width)) / sponge_width)^2
    z_mask = max(0, (z + top_sponge) / top_sponge)^2
    return max(x_mask, z_mask)
end

@inline u_damping(x, y, z, t, u) = - sponge_mask(x, z) * (u - U0) / 600
@inline u_damping(x, z, t, u)    = - sponge_mask(x, z) * (u - U0) / 600
@inline u_damping(x, y, z, t)    = - sponge_mask(x, z) * 0.0
@inline u_damping(x, z, t)       = - sponge_mask(x, z) * 0.0

@inline v_damping(x, y, z, t, v) = - sponge_mask(x, z) * v / 600
@inline v_damping(x, z, t, v)    = - sponge_mask(x, z) * v / 600
@inline v_damping(x, y, z, t)    = - sponge_mask(x, z) * 0.0
@inline v_damping(x, z, t)       = - sponge_mask(x, z) * 0.0

@inline w_damping(x, y, z, t, w) = - sponge_mask(x, z) * w / 600
@inline w_damping(x, z, t, w)    = - sponge_mask(x, z) * w / 600
@inline w_damping(x, y, z, t)    = - sponge_mask(x, z) * 0.0
@inline w_damping(x, z, t)       = - sponge_mask(x, z) * 0.0

u_forcing = Forcing(u_damping, field_dependencies = :u)
v_forcing = Forcing(v_damping, field_dependencies = :v)
w_forcing = Forcing(w_damping, field_dependencies = :w)

# 5. Model Setup
pressure_solver = ConjugateGradientPoissonSolver(grid)

model = NonhydrostaticModel(
    grid,
    pressure_solver = pressure_solver,
    advection = Centered(),
    timestepper = :RungeKutta3,
    coriolis = f0 == 0 ? nothing : FPlane(f = f0),   # CHANGED: allow f0 = 0
    buoyancy = BuoyancyTracer(),
    tracers = :b,
    background_fields = (b = b_field,),
    forcing = (u = u_forcing, v = v_forcing, w = w_forcing)
)

set!(model, u = U0)

# 6. Simulation & Output
simulation = Simulation(model, Δt = 20seconds, stop_time = 18hours)   # CHANGED: Δt was 5 s

progress(sim) = @printf("Time: %s, Iteration: %d\n", prettytime(sim), iteration(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(200))

simulation.output_writers[:fields] = JLD2Writer(
    model, model.velocities,
    filename = "internal_waves_corrected.jld2",   # CHANGED name
    schedule = TimeInterval(10minutes),
    overwrite_files = true
)

run!(simulation)

# 7. Visualization: Vertical Velocity Wave Field
println("Generating wave animation...")

u_ts = FieldTimeSeries("internal_waves_corrected.jld2", "u")   # CHANGED name
w_ts = FieldTimeSeries("internal_waves_corrected.jld2", "w")   # CHANGED name

times = w_ts.times
n_frames = length(times)

xc, yc, zc = nodes(w_ts)

n = Observable(1)
w_frame = @lift interior(w_ts[$n], :, 1, :)

fig1 = Figure(size = (900, 450))
ax1 = Axis(fig1[1, 1],
           title = @lift("Vertical Velocity w [m/s] at t = " * prettytime(times[$n])),
           xlabel = "x [km]", ylabel = "z [m]")

hm = heatmap!(ax1, xc ./ 1000, zc, w_frame, colorrange = (-0.005, 0.005), colormap = :balance)
Colorbar(fig1[1, 2], hm, label = "w (m/s)")

x_topo = range(-Lx/2, Lx/2, length=500)
lines!(ax1, x_topo ./ 1000, bottom_topography.(x_topo), color = :black, linewidth = 2)

record(fig1, "internal_wave_generation_corrected.gif", 1:n_frames, framerate = 12) do i   # CHANGED name
    n[] = i
end

println("Done! Animation saved as internal_wave_generation_corrected.gif")

# 8. Quantitative Diagnostic: Wave Drag & Momentum Flux vs. Linear Theory
println("Computing momentum flux diagnostic...")

D_lin = N * U0 * h0^2  # Linear theory hydrostatic drag per unit density
zw = zc
k_idx = argmin(abs.(zw .+ 600))  # Mid-depth level below the top sponge
dx = Lx / Nx

flux = map(1:length(w_ts.times)) do i
    u_data = interior(u_ts[i], :, 1, k_idx)
    u_center = (u_data .+ circshift(u_data, -1)) ./ 2
    w_data = interior(w_ts[i], :, 1, k_idx)
    return sum((u_center .- U0) .* w_data) * dx
end

fig2 = Figure(size = (650, 400))
ax2 = Axis(fig2[1, 1],
           title = "Vertical Momentum Flux ∫u'w' dx at z = -600 m",
           xlabel = "Time [hours]", ylabel = "Flux [m³/s²]")

lines!(ax2, times ./ 3600, flux, label = "Simulation (Oceananigans)", color = :blue, linewidth = 2)
hlines!(ax2, [-D_lin], linestyle = :dash, color = :red, label = "Linear Theory (-N U₀ h₀²)")
axislegend(ax2, position = :rt)   # CHANGED: was :rb, which covered the curve

save("flux_vs_theory_corrected.png", fig2)   # CHANGED name
println("Done! Quantitative plot saved as flux_vs_theory_corrected.png")
