# =====================================================================
# replot_demo_flux.jl  --  redraw the flux plot of internal_wave_demo_corrected.jl
# from its saved output, without re-running the simulation.
#
# Usage: put this file in the same folder as internal_waves_corrected.jld2, then
#   julia --project replot_demo_flux.jl
#
# It uses exactly the diagnostic of the demo script (flux at z = -600 m, integrated over x),
# draws the legend at the top right, writes flux_vs_theory_corrected.png, and prints the mean
# drag over the last 2 h as a ratio to N U h0^2 together with its drift over that window.
# =====================================================================
cd(@__DIR__)
using Oceananigans, CairoMakie, Statistics

u_ts = FieldTimeSeries("internal_waves_corrected.jld2", "u")
w_ts = FieldTimeSeries("internal_waves_corrected.jld2", "w")
times = w_ts.times
xc, yc, zc = nodes(w_ts)
k_idx = argmin(abs.(zc .+ 600))

dx = 20e3 / 512
U0 = 0.2
N = sqrt(4e-6)
h0 = 20.0
D_lin = N * U0 * h0^2

flux = map(1:length(times)) do i
    u_data = interior(u_ts[i], :, 1, k_idx)
    u_center = (u_data .+ circshift(u_data, -1)) ./ 2
    w_data = interior(w_ts[i], :, 1, k_idx)
    sum((u_center .- U0) .* w_data) * dx
end

fig2 = Figure(size = (650, 400))
ax2 = Axis(fig2[1, 1], title = "Vertical Momentum Flux ∫u'w' dx at z = -600 m",
           xlabel = "Time [hours]", ylabel = "Flux [m³/s²]")
lines!(ax2, times ./ 3600, flux, label = "Simulation (Oceananigans)", color = :blue, linewidth = 2)
hlines!(ax2, [-D_lin], linestyle = :dash, color = :red, label = "Linear Theory (-N U₀ h₀²)")
axislegend(ax2, position = :rt)
save("flux_vs_theory_corrected.png", fig2)
println("Saved flux_vs_theory_corrected.png")

win = findall(t -> t >= times[end] - 2 * 3600, times)
println("mean drag over the last 2 h / (N U h0^2) = ", round(-mean(flux[win]) / D_lin, digits = 3), "   (linear theory: 0.980)")
println("drift over that window = ", round(100 * (flux[win[end]] - flux[win[1]]) / mean(flux[win]), digits = 1), " %")
