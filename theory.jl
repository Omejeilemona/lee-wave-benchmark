# =====================================================================
# theory.jl  --  the "answer key": linear lee-wave theory for the ridge
# =====================================================================
# Needs only FFTW (no Oceananigans). Run:   julia --project theory.jl
#
# SOURCES
#   Linear steady lee-wave theory: Bell (1975), J. Geophys. Res. 80, 320-327.
#   Formulation used here (equation numbers refer to the arXiv preprint,
#   please check them against the published version before quoting):
#     Baker & Mashayek (2021), "Surface reflection of bottom generated
#     oceanic lee waves", J. Fluid Mech. (arXiv:2103.03779).
#       eq 2.18  vertical wavenumber m(k)
#       eq 2.19  waves propagate only if |f| < |U k| < |N|
#       eq 2.21-2.23  bottom condition w = U dh/dx, solution via Fourier transform
#       eq 3.1   radiating solution  zeta_hat = exp(i m z)
#       eq 2.8, 2.10  give v and b from u and w (for uniform U, no viscosity)
#       eq 2.46  Eliassen-Palm flux  F = <u w> - f <v b> / N^2
#       eq 2.48  energy flux  <p w> = -rho0 * U * F
#   This is a fresh implementation of those equations (not their code).
#
# WHAT IT COMPUTES
#   * the wave drag per unit density,  D = -F  [m^3/s^2], with and without
#     rotation, integrated over x;
#   * the linear vertical-velocity field w(x, z) to compare with the
#     simulation.
#
# Sign convention: for U > 0 upward-radiating waves have <u w> < 0, so F < 0
# and the drag D = -F is positive.
#
# Why F and not just <u w>: with rotation, <u w> is not conserved with height
# (Bretherton 1969); F is (Baker & Mashayek 2021, sec 2.8).
# =====================================================================

using FFTW
using Printf
include(joinpath(@__DIR__, "params.jl"))

# ---- vertical wavenumber m(k), eq 2.18 --------------------------------
# m^2 = k^2 (N^2 - alpha U^2 k^2) / (U^2 k^2 - f^2)
# Radiating (wave-like) if |f| < |U k| < |N| (nonhydrostatic: alpha = 1);
# otherwise the disturbance decays with height (m imaginary, decaying root).
# For radiating waves m takes the sign of U*k so energy goes upward.
function vertical_wavenumber(k, N, U, f, α)
    denom = U^2 * k^2 - f^2
    denom == 0 && return 0.0 + 0.0im        # degenerate single wavenumber, ignored
    m² = k^2 * (N^2 - α * U^2 * k^2) / denom
    if denom > 0 && N^2 > α * U^2 * k^2
        return sign(U * k) * sqrt(m²) + 0.0im
    else
        return 0.0 + im * sqrt(abs(m²))
    end
end

# ---- ridge and its Fourier transform ----------------------------------
# A periodic box much wider than the ridge (400 km) so the ridge behaves
# as an isolated one.
function ridge_spectrum(p; L = 400e3, n = 2^15)
    dx = L / n
    x  = ((0:n-1) .- n / 2) .* dx                 # x = 0 at the centre
    h  = @. p.h0 * exp(-x^2 / (2 * p.σ^2))
    ĥ  = fft(ifftshift(h)) .* dx                   # continuous-transform approximation
    k  = 2π .* fftfreq(n, 1 / dx)                  # rad/m
    k  = [ki == 0 ? 1e-30 : ki for ki in k]        # avoid 0/0 at k = 0
    return (; x, h, ĥ, k, dx, L, n)
end

# ---- wave fields in wavenumber space at height z above the ridge base --
# psi_hat = U h_hat exp(i m z)            (eq 2.22, 2.23, 3.1)
# w = psi_x,  u = -psi_z                   (definitions below eq 2.11)
# v_hat = i f u_hat/(k U)  from  U v_x + f u = 0   (eq 2.8)
# b_hat = i N^2 w_hat/(k U) from  U b_x + N^2 w = 0 (eq 2.10, U_z = 0)
function wave_fields_hat(p, S, z)
    m  = vertical_wavenumber.(S.k, p.N, p.U, p.f, p.α)
    ψ̂  = @. p.U * S.ĥ * exp(im * m * z)
    ŵ  = @. im * S.k * ψ̂
    û  = @. -im * m * ψ̂
    v̂  = @. im * p.f * û / (S.k * p.U)
    b̂  = @. im * p.N^2 * ŵ / (S.k * p.U)
    return (; û, v̂, ŵ, b̂)
end

# Parseval: integral of a*c over x  =  sum(a_hat * conj(c_hat)) / L
parseval(â, ĉ, S) = real(sum(â .* conj.(ĉ))) / S.L

"""
    linear_fluxes(p; z)

Integrated (over x) wave fluxes at height z above the ridge base, per unit
density [m^3/s^2]. Returns uw, vb, F = uw - f*vb/N^2 and drag = -F.
"""
function linear_fluxes(p; z)
    S  = ridge_spectrum(p)
    q  = wave_fields_hat(p, S, z)
    uw = parseval(q.û, q.ŵ, S)
    vb = parseval(q.v̂, q.b̂, S)
    F  = uw - p.f * vb / p.N^2
    return (; uw, vb, F, drag = -F)
end

"""
    linear_w(p; z)

Linear vertical-velocity field w(x) at height z above the ridge base.
Returns (x, w).
"""
function linear_w(p; z)
    S = ridge_spectrum(p)
    q = wave_fields_hat(p, S, z)
    w = fftshift(real(ifft(q.ŵ)) ./ S.dx)
    return S.x, w
end

# Hydrostatic, non-rotating reference value for a Gaussian ridge:
# drag = N*U*h0^2  (the number in your original script).
reference_drag(p) = p.N * p.U * p.h0^2

# ---- self-test: run `julia --project theory.jl` ------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    for f in (0.0, f_arctic)
        p = benchmark_params(; f)
        describe(p)
        ref = reference_drag(p)
        r = linear_fluxes(p; z = 100.0)
        @printf("reference N*U*h0^2           = %.4f m^3/s^2\n", ref)
        @printf("-integral(u'w')/ref          = %.3f   (momentum flux alone)\n", -r.uw / ref)
        @printf("drag = -F / ref              = %.3f   (conserved flux; compare the simulation to this)\n", r.drag / ref)

        # sanity check: bottom boundary condition w(x, 0) = U dh/dx
        x, w0 = linear_w(p; z = 0.0)
        dhdx  = @. -x / p.σ^2 * p.h0 * exp(-x^2 / (2 * p.σ^2))
        err   = maximum(abs.(w0 .- p.U .* dhdx)) / maximum(abs.(p.U .* dhdx))
        @printf("check w(x,0) = U dh/dx: relative error = %.1e   (should be ~1e-4 or smaller)\n\n", err)
    end
end

# EXPECTED OUTPUT (from a Python port of this logic, N = 2e-3, U = 0.2, h0 = 20, sigma = 500):
#   f = 0       : -integral(u'w')/ref = 0.980 ; drag/ref = 0.980
#   f = 1.4e-4  : -integral(u'w')/ref = 1.070 ; drag/ref = 0.748
#   w(x,0) check ~ 2e-4
# If your Julia numbers differ, tell me before running the simulation.
