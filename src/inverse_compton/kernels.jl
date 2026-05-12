using NumericalIntegration


"""
    u_CMB(z)

CMB photon energy density at redshift `z` [erg cm⁻³].
"""
u_CMB(z::Real) = (4σ_SB / cL) * (T_CMB0 * (1.0 + z))^4


"""
    ε_CMB(z)

Mean CMB photon energy at redshift `z` [erg], using ε₀ = 2.7 k_B T.
"""
ε_CMB(z::Real) = ic_mean_photon_factor * k_B * T_CMB0 * (1.0 + z)


"""
    ic_kernel(q)

Thomson IC kernel from Blumenthal & Gould (1970), valid for 0 < q ≤ 1,
where q = E_X / (4γ² ε₀).
"""
@inline function ic_kernel(q::Real)
    (q <= 0.0 || q > 1.0) && return 0.0
    return 2q * log(q) + (1 + 2q) * (1 - q)
end


"""
    n_CMB_planck(ε, T)

Planck differential photon number density at energy `ε` [erg] and
temperature `T` [K], in [photons cm⁻³ erg⁻¹]:

    n(ε, T) = (8π / (h c)³) · ε² / (exp(ε / k_B T) − 1)

Returns 0 for ε ≤ 0. Uses `expm1` to avoid catastrophic cancellation in
the Rayleigh-Jeans regime.
"""
@inline function n_CMB_planck(ε::Real, T::Real)
    ε <= 0.0 && return 0.0
    x = ε / (k_B * T)
    # expm1(x) = exp(x) - 1, accurate near x → 0
    return ic_planck_prefac * ε^2 / expm1(x)
end


"""
    ic_planck_photon_integral(γ², E_X, T; N=N_planck_default)

Inner photon-bath integral for the Planck-CMB IC source:

    I_ph(γ, E_X) = ∫_{ε_min}^{ε_max} dε  n_CMB(ε, T)/ε  ·  g(E_X/(4γ² ε))

where `ε_min = E_X / (4γ²)` enforces q ≤ 1 in the Thomson kernel `g`.
Upper limit `ε_max = max(30 k_B T, 100 ε_min)` covers both the Wien tail
when ε_min ≲ k_B T and the exponentially-suppressed regime when
ε_min ≫ k_B T. Trapezoidal on a log-spaced ε grid.

Returns [cm⁻³ erg⁻¹] — same units as `u_CMB / ε₀²` in the monochromatic
branch, so it slots into `ic_source_per_bin` with the existing
`ic_prefac_c` prefactor unchanged.
"""
function ic_planck_photon_integral(γ2::Real, E_X::Real, T::Real;
                                    N::Int=N_planck_default)
    ε_min = E_X / (4.0 * γ2)
    ε_max = max(30.0 * k_B * T, 100.0 * ε_min)
    (ε_max <= ε_min) && return 0.0

    ε_grid = 10.0 .^ LinRange(log10(ε_min), log10(ε_max), N)
    integrand = Vector{Float64}(undef, N)
    @inbounds for i = 1:N
        ε = ε_grid[i]
        q = E_X / (4.0 * γ2 * ε)
        integrand[i] = n_CMB_planck(ε, T) / ε * ic_kernel(q)
    end

    return integrate(ε_grid, integrand, Trapezoidal())
end
