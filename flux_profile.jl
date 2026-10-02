# =====================================================================
# flux_profile.jl  --  wave flux versus HEIGHT in the simulation (diagnostic)
# =====================================================================
# Usage (inside this folder, after a run has finished):
#   run(`julia --project flux_profile.jl lee_f0_nz512.jld2 0.0`)
#   run(`julia --project flux_profile.jl lee_f0.jld2 0.0`)
#
# It averages the drag  -F = -(u'w' - f v'b'/N^2), integrated over x outside the
# sponges, over t = 4-6 h, at several heights above the ridge crest.
#
# For f = 0 (and f = 1.4e-4 up to 300 m) it prints, next to the simulation, my linear PREDICTION for the same
# quantity: an impulsive start of flow over the ridge (inviscid, nonhydrostatic,
# f = 0, N = 2e-3, U = 0.2, h0 = 20, sigma = 500, a tall rigid-lid domain),
# averaged over t = 4-6 h. This is my own calculation (a small Python script),
# not a published result. It uses the same equations as theory.jl, but marched in
# time from the start-up, so it includes the fact that the waves have not yet
# reached the upper levels.
#
# How to read it:
#   * ratio simulation/linear about the same at all heights  -> the ridge is forcing
#     waves with the wrong amplitude (geometry, boundary treatment or nonlinearity).
#   * ratio changes strongly with height -> the problem is in how the waves propagate
#     (numerical dispersion/dissipation, or a different effective U or N).
# =====================================================================

using Oceananigans
using Printf
using Statistics
include(joinpath(@__DIR__, "params.jl"))

file  = ARGS[1]
f_val = parse(Float64, ARGS[2])

u_ts = FieldTimeSeries(file, "u")
v_ts = FieldTimeSeries(file, "v")
w_ts = FieldTimeSeries(file, "w")
b_ts = FieldTimeSeries(file, "b")
times = u_ts.times

get2d(ts, n) = Array(interior(ts[n], :, 1, :))
Nx_file, Nz_file = size(get2d(u_ts, 1))
m_h0   = match(r"_h(\d+)", basename(file))          # ridge height is encoded in the file name, e.g. _h40
h0_file = m_h0 === nothing ? 20.0 : parse(Float64, m_h0.captures[1])
p = benchmark_params(; f = f_val, Nx = Nx_file, Nz = Nz_file, h0 = h0_file)

dx = p.Lx / p.Nx
dz = p.Lz / p.Nz
xc = [-p.Lx / 2 + (i - 0.5) * dx for i in 1:p.Nx]
zc = [-p.Lz + (k - 0.5) * dz for k in 1:p.Nz]
crest  = -p.Lz + p.h0
inside = findall(x -> abs(x) < p.Lx / 2 - p.sponge_width, xc)

center_u(A) = (A .+ circshift(A, (-1, 0))) ./ 2
center_w(A) = (A[:, 1:end-1] .+ A[:, 2:end]) ./ 2

heights = [25, 50, 75, 100, 150, 200, 300, 400, 500]       # m above the crest
kidx    = [argmin(abs.(zc .- (crest + h))) for h in heights]
win     = findall(t -> 4 * 3600 - 1 <= t <= 6 * 3600 + 1, times)
isempty(win) && error("No snapshots between 4 and 6 h in this file.")

drag = zeros(length(heights))
for n in win
    u = center_u(get2d(u_ts, n))
    w = center_w(get2d(w_ts, n))
    v = get2d(v_ts, n)
    b = get2d(b_ts, n)
    for (m, k) in enumerate(kidx)
        uw = sum((u[inside, k] .- p.U) .* w[inside, k]) * dx
        vb = sum(v[inside, k] .* b[inside, k]) * dx
        drag[m] += -(uw - p.f * vb / p.N^2) / length(win)
    end
end
ref = p.N * p.U * p.h0^2      # N*U*h0^2, the same reference as in theory.jl

# linear impulsive-start prediction, mean over 4-6 h (f = 0 only; see header)
linear = [1.032, 1.018, 0.997, 0.972, 0.909, 0.831, 0.674, 0.527, 0.404]
# same model with rotation f = 1.4e-4 (first 7 heights only: 25-300 m), mean over 4-6 h
linear_rot = [0.686, 0.694, 0.689, 0.677, 0.640, 0.590, 0.482]

println("\n=========== flux versus height (mean over t = 4-6 h, f = $(p.f), Nz = $(p.Nz)) ===========")
if p.f == 0
    println("height above crest [m]   simulation/ref   linear/ref   ratio sim/linear")
    for (m, h) in enumerate(heights)
        @printf("%8d                  %7.3f        %7.3f      %6.2f\n", h, drag[m] / ref, linear[m], drag[m] / ref / linear[m])
    end
else
    println("height above crest [m]   simulation/ref   linear/ref   ratio sim/linear   (rotating linear values stored for f = 1.4e-4 only)")
    for (m, h) in enumerate(heights)
        if m <= length(linear_rot)
            @printf("%8d                  %7.3f        %7.3f      %6.2f\n", h, drag[m] / ref, linear_rot[m], drag[m] / ref / linear_rot[m])
        else
            @printf("%8d                  %7.3f            -           -\n", h, drag[m] / ref)
        end
    end
end
println("=====================================================================================\n")
