# =====================================================================
# compare.jl  --  simulation vs linear theory
# =====================================================================
# Usage from the julia> prompt, inside this folder:
#
#   run(`julia --project compare.jl lee_f0.jld2 0.0`)
#   run(`julia --project compare.jl lee_frot.jld2 1.4e-4`)
#
# (first argument = output file from run_case.jl, second = the SAME f used for the run)
#
# What it does
#   1. Computes the wave drag from the simulation at 50-100 m above the ridge crest,
#      integrated over the interior of the domain (outside the sponges), and averaged
#      over the last 2 hours of output.
#   2. Compares it with linear theory (theory.jl), both for F = u'w' - f v'b'/N^2
#      (the conserved flux) and for u'w' alone.
#   3. Saves two figures: fig_w_<tag>.png and fig_flux_<tag>.png.
# =====================================================================

using Oceananigans
using CairoMakie
using Printf
using Statistics
include(joinpath(@__DIR__, "theory.jl"))

file  = ARGS[1]
f_val = parse(Float64, ARGS[2])
tag   = replace(basename(file), ".jld2" => "")

u_ts = FieldTimeSeries(file, "u")
v_ts = FieldTimeSeries(file, "v")
w_ts = FieldTimeSeries(file, "w")
b_ts = FieldTimeSeries(file, "b")
times = u_ts.times

# grid size is read from the file, so runs at different resolutions compare correctly
Nx_file, Nz_file = size(Array(interior(u_ts[1], :, 1, :)))
m_h0   = match(r"_h(\d+)", basename(file))          # ridge height is encoded in the file name, e.g. _h40
h0_file = m_h0 === nothing ? 20.0 : parse(Float64, m_h0.captures[1])
p = benchmark_params(; f = f_val, Nx = Nx_file, Nz = Nz_file, h0 = h0_file)

# ---- grid coordinates (uniform grid, so computed directly) ---------------
dx = p.Lx / p.Nx
dz = p.Lz / p.Nz
xc = [-p.Lx / 2 + (i - 0.5) * dx for i in 1:p.Nx]
zc = [-p.Lz + (k - 0.5) * dz for k in 1:p.Nz]       # cell centres
zf = [-p.Lz + (k - 1) * dz for k in 1:p.Nz+1]       # cell faces (w lives here)

crest  = -p.Lz + p.h0
ks     = findall(z -> crest + p.z_above_crest[1] <= z <= crest + p.z_above_crest[2], zc)
inside = findall(x -> abs(x) < p.Lx / 2 - p.sponge_width, xc)   # keep away from the x-sponges
isempty(ks) && error("No grid levels found in the measurement band.")

# ---- staggered-grid helpers (move everything to cell centres) -----------
center_u(A) = (A .+ circshift(A, (-1, 0))) ./ 2          # x-faces -> x-centres (periodic in x)
center_w(A) = (A[:, 1:end-1] .+ A[:, 2:end]) ./ 2        # z-faces -> z-centres
get2d(ts, n) = Array(interior(ts[n], :, 1, :))

# ---- fluxes from the simulation, integrated over x ----------------------
function snapshot_fluxes(n)
    u = center_u(get2d(u_ts, n))
    w = center_w(get2d(w_ts, n))
    v = get2d(v_ts, n)
    b = get2d(b_ts, n)
    uw = mean([sum((u[inside, k] .- p.U) .* w[inside, k]) * dx for k in ks])
    vb = mean([sum(v[inside, k] .* b[inside, k]) * dx for k in ks])
    return uw, vb
end

nt  = length(times)
uw  = zeros(nt)
vb  = zeros(nt)
for n in 1:nt
    uw[n], vb[n] = snapshot_fluxes(n)
end
F    = uw .- p.f .* vb ./ p.N^2
drag = .-F

win = findall(t -> t >= times[end] - p.avg_window, times)
ref = reference_drag(p)
th  = linear_fluxes(p; z = 100.0)

sim_drag   = mean(drag[win])
sim_uw     = mean(.-uw[win])

println("\n================ comparison (f = $(p.f)) ================")
@printf("snapshots used for the average: %d  (t = %.2f to %.2f h)\n", length(win), times[win[1]] / 3600, times[win[end]] / 3600)
@printf("reference N*U*h0^2                     = %.4f m^3/s^2\n", ref)
@printf("                          theory/ref   simulation/ref   difference\n")
@printf("drag = -F (conserved)     %8.3f     %8.3f        %+6.1f %%\n", th.drag / ref, sim_drag / ref, 100 * (sim_drag - th.drag) / th.drag)
@printf("-integral u'w' only       %8.3f     %8.3f        %+6.1f %%\n", -th.uw / ref, sim_uw / ref, 100 * (sim_uw + th.uw) / (-th.uw))
@printf("drift of drag over the averaging window: %.1f %% of its mean\n", 100 * (drag[win[end]] - drag[win[1]]) / sim_drag)
println("=======================================================\n")

# ---- Figure 1: vertical velocity -----------------------------------------
k0 = argmin(abs.(zf .- (crest + 75)))                    # w level 75 m above the crest
w_last = get2d(w_ts, nt)
w_sim_line = mean([center_w(get2d(w_ts, n))[:, argmin(abs.(zc .- (crest + 75)))] for n in win])
x_th, w_th = linear_w(p; z = p.h0 + 75.0)               # height above ridge base
mask_th = findall(x -> abs(x) < p.Lx / 2 - p.sponge_width, x_th)
wlim = 1.5 * maximum(abs, w_th)

fig1 = Figure(size = (1000, 750))
ax1 = Axis(fig1[1, 1], title = "Simulated w at t = $(round(times[end] / 3600, digits = 1)) h  (f = $(p.f))",
           xlabel = "x [km]", ylabel = "z [m]")
hm = heatmap!(ax1, xc ./ 1000, zf, w_last, colorrange = (-wlim, wlim), colormap = :balance)
Colorbar(fig1[1, 2], hm, label = "w [m/s]")
vlines!(ax1, [-(p.Lx / 2 - p.sponge_width) / 1000, (p.Lx / 2 - p.sponge_width) / 1000], color = :gray, linestyle = :dash)
hlines!(ax1, [-p.top_sponge], color = :gray, linestyle = :dash)

ax2 = Axis(fig1[2, 1], title = "w at 75 m above the crest: simulation (time-mean) vs linear theory",
           xlabel = "x [km]", ylabel = "w [m/s]")
lines!(ax2, xc[inside] ./ 1000, w_sim_line[inside], label = "simulation", linewidth = 2)
lines!(ax2, x_th[mask_th] ./ 1000, w_th[mask_th], label = "linear theory", linewidth = 2, linestyle = :dash, color = :red)
axislegend(ax2, position = :rt)
save("fig_w_$(tag).png", fig1)

# ---- Figure 2: drag vs time ----------------------------------------------
fig2 = Figure(size = (800, 450))
ax3 = Axis(fig2[1, 1], title = "Wave drag at 50-100 m above the crest  (f = $(p.f))",
           xlabel = "time [h]", ylabel = "drag / (N U h0^2)")
lines!(ax3, times ./ 3600, drag ./ ref, linewidth = 2, label = "simulation: -F")
lines!(ax3, times ./ 3600, (.-uw) ./ ref, linewidth = 1, linestyle = :dot, label = "simulation: -u'w' only")
hlines!(ax3, [th.drag / ref], color = :red, linestyle = :dash, label = "theory: -F")
p.f != 0 && hlines!(ax3, [-th.uw / ref], color = :gray, linestyle = :dash, label = "theory: -u'w' only")
axislegend(ax3, position = :rb)
save("fig_flux_$(tag).png", fig2)

println("Saved fig_w_$(tag).png and fig_flux_$(tag).png")
