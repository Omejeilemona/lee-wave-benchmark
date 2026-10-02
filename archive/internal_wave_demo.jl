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
h0 = 20meters          # Linear regime: N*h0/U ~ 0.3 << 1
sigma = 500meters      # Narrow ridge: Uk ~ 4e-4 s⁻¹ > f0

bottom_topography(x) = -Lz + h0 * exp(-x^2 / (2 * sigma^2))

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_topography))

# 3. Arctic Stratification, Coriolis & Flow Parameters
N2 = 1e-4              # N ~ 0.01 s⁻¹
N = sqrt(N2)
U0 = 0.2               # Background inflow velocity (m/s)
f0 = 1.4e-4            # Arctic Coriolis parameter (~80°N)

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
    coriolis = FPlane(f = f0),
    buoyancy = BuoyancyTracer(),
    tracers = :b,
    background_fields = (b = b_field,),
    forcing = (u = u_forcing, v = v_forcing, w = w_forcing)
)

set!(model, u = U0)

# 6. Simulation & Output
simulation = Simulation(model, Δt = 5seconds, stop_time = 18hours)

progress(sim) = @printf("Time: %s, Iteration: %d\n", prettytime(sim), iteration(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(200))

simulation.output_writers[:fields] = JLD2Writer(
    model, model.velocities,
    filename = "internal_waves.jld2",
    schedule = TimeInterval(10minutes),
    overwrite_files = true
)

run!(simulation)

# 7. Visualization: Vertical Velocity Wave Field
println("Generating wave animation...")

u_ts = FieldTimeSeries("internal_waves.jld2", "u")
w_ts = FieldTimeSeries("internal_waves.jld2", "w")

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

record(fig1, "internal_wave_generation.gif", 1:n_frames, framerate = 12) do i
    n[] = i
end

println("Done! Animation saved as internal_wave_generation.gif")

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
axislegend(ax2, position = :rb)

save("flux_vs_theory.png", fig2)
println("Done! Quantitative plot saved as flux_vs_theory.png")
