using NumericalIntegration


"""
IC contribution from one momentum bin using 3-point Simpson's rule.

`photon_factor(γ², E_X)` returns the photon-bath factor in [cm⁻³ erg⁻¹]:
- monochromatic: `(γ², E_X) -> ic_kernel(E_X / (4 γ² ε₀))` (the `u/ε₀²`
  normalization is folded into the outer prefactor),
- Planck: `(γ², E_X) -> ic_planck_photon_integral(γ², E_X, T; N)`
  (carries its own normalization, outer prefactor is just `ic_prefac_c`).
"""
function ic_source_per_bin(f_p_start::Real, p_start::Real, p_end::Real,
                            q_slope::Real, E_X::Real,
                            photon_factor::F) where {F}

    γ2 = p_start^2 + 1.0
    F_start = p_start^2 * f_p_start / γ2 * photon_factor(γ2, E_X)

    p_mid = find_log_mid(p_start, p_end)
    f_mid = interpolate_spectrum(p_mid, f_p_start, p_start, q_slope)
    γ2 = p_mid^2 + 1.0
    F_mid = p_mid^2 * f_mid / γ2 * photon_factor(γ2, E_X)

    f_end = interpolate_spectrum(p_end, f_p_start, p_start, q_slope)
    γ2 = p_end^2 + 1.0
    F_end = p_end^2 * f_end / γ2 * photon_factor(γ2, E_X)

    return (p_end - p_start) / 6.0 * (F_start + 4F_mid + F_end)
end


"""
    ic_source_cmb(f_p, q, cut, bounds, E_X; z, reduce_spectrum, N_subcycle,
                  planck_cmb, N_planck)

IC source function [photons s⁻¹ erg⁻¹ cm⁻³] at X-ray photon energy `E_X` [erg]
for CR electrons upscattering CMB photons at redshift `z`.

# Arguments
- `f_p`    : spectral norms at bin starts
- `q`      : power-law slopes per bin
- `cut`    : cutoff momentum [m_e c]
- `bounds` : N+1 bin boundaries [m_e c]
- `E_X`    : photon energy [erg]

# Keyword Arguments
- `z::Real=0.0`                : redshift
- `reduce_spectrum::Bool=true` : if false, return (bin_centers, per-bin values)
- `N_subcycle::Int=10`         : subcycle bins for accuracy (1 = single Simpson step)
- `planck_cmb::Bool=false`     : if true, integrate over the full Planck CMB
                                 spectrum instead of using the monochromatic
                                 ε₀ = 2.7 k_B T approximation
- `N_planck::Int=N_planck_default` : number of points for the inner ε-integration
                                     when `planck_cmb=true` (ignored otherwise)
"""
function ic_source_cmb(f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        bounds::Vector{<:Real},
                        E_X::Real;
                        z::Real=0.0,
                        reduce_spectrum::Bool=true,
                        N_subcycle::Int=10,
                        planck_cmb::Bool=false,
                        N_planck::Int=N_planck_default)

    if iszero(sum(f_p))
        return 0.0
    end

    Nbins = length(f_p)

    if planck_cmb
        T = T_CMB0 * (1.0 + z)
        prefac = ic_prefac_c
        photon_factor = (γ2, EX) -> ic_planck_photon_integral(γ2, EX, T; N=N_planck)
        # In Planck mode, the effective max photon energy is the Wien cap that
        # bounds the inner integral. Below this, even the hardest CMB photon
        # in the integration range cannot upscatter to E_X.
        ε_eff_max = 30.0 * k_B * T
        p_min_ic = sqrt(max(E_X / (4.0 * ε_eff_max) - 1.0, 0.0))
    else
        ε0     = ε_CMB(z)
        prefac = ic_prefac_c * u_CMB(z) / ε0^2
        photon_factor = (γ2, EX) -> ic_kernel(EX / (4.0 * γ2 * ε0))
        # minimum momentum to upscatter to E_X: E_X = 4γ²ε₀ → p_min = √(E_X/4ε₀ - 1)
        p_min_ic = sqrt(max(E_X / (4.0 * ε0) - 1.0, 0.0))
    end

    q_ic = Vector{Float64}(undef, Nbins)

    if !reduce_spectrum
        bin_centers = Vector{Float64}(undef, Nbins)
    end

    @inbounds for i = 1:Nbins

        p_start = bounds[i]

        if p_start > cut
            q_ic[i:end] .= 0.0
            break
        end

        p_end = bounds[i+1]
        if p_end > cut
            p_end = cut
        end

        if !reduce_spectrum
            bin_centers[i] = find_log_mid(p_start, p_end)
        end

        # all electrons in this bin too slow to reach E_X
        if p_end < p_min_ic
            q_ic[i] = 0.0
            continue
        end

        f_p_start = f_p[i]

        # threshold falls within bin: advance start to p_min_ic
        if p_start < p_min_ic
            f_p_start = interpolate_spectrum(p_min_ic, f_p_start, p_start, q[i])
            p_start = p_min_ic
        end

        if isnan(f_p_start)
            q_ic[i] = 0.0
            continue
        end

        if N_subcycle > 1
            q_ic[i] = 0.0
            p_sub = 10.0 .^ LinRange(log10(p_start), log10(p_end), N_subcycle + 1)
            for j = 1:N_subcycle
                f_sub = interpolate_spectrum(p_sub[j], f_p_start, p_start, q[i])
                q_ic[i] += ic_source_per_bin(f_sub, p_sub[j], p_sub[j+1], q[i], E_X, photon_factor)
            end
        else
            q_ic[i] = ic_source_per_bin(f_p_start, p_start, p_end, q[i], E_X, photon_factor)
        end

    end

    if reduce_spectrum
        return prefac * sum(q_ic)
    else
        return bin_centers, prefac .* q_ic
    end
end


"""
    ic_source_cmb(f_p, q, cut, par, E_X; z, reduce_spectrum, N_subcycle,
                  planck_cmb, N_planck)
"""
function ic_source_cmb(f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        par::CRMomentumDistributionConfig,
                        E_X::Real;
                        z::Real=0.0,
                        reduce_spectrum::Bool=true,
                        N_subcycle::Int=10,
                        planck_cmb::Bool=false,
                        N_planck::Int=N_planck_default)

    bounds = momentum_bin_boundaries(par)
    return ic_source_cmb(f_p, q, cut, bounds, E_X;
                         z, reduce_spectrum, N_subcycle, planck_cmb, N_planck)
end


"""
    ic_source_cmb(CRe, par, E_X; z, reduce_spectrum, N_subcycle,
                  planck_cmb, N_planck)
"""
function ic_source_cmb(CRe::CRMomentumDistribution,
                        par::CRMomentumDistributionConfig,
                        E_X::Real;
                        z::Real=0.0,
                        reduce_spectrum::Bool=true,
                        N_subcycle::Int=10,
                        planck_cmb::Bool=false,
                        N_planck::Int=N_planck_default)

    f_p, q, cut = convert(CRe)
    bounds = momentum_bin_boundaries(par)
    return ic_source_cmb(f_p, q, cut, bounds, E_X;
                         z, reduce_spectrum, N_subcycle, planck_cmb, N_planck)
end


"""
    ic_emissivity_cmb(f_p, q, cut, bounds, E_X; z, N_subcycle,
                      planck_cmb, N_planck)

IC emissivity `E_X × source` [erg s⁻¹ erg⁻¹ cm⁻³] at photon energy `E_X` [erg].
See `ic_source_cmb` for the argument list.
"""
function ic_emissivity_cmb(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            E_X::Real;
                            z::Real=0.0,
                            N_subcycle::Int=10,
                            planck_cmb::Bool=false,
                            N_planck::Int=N_planck_default)

    return E_X * ic_source_cmb(f_p, q, cut, bounds, E_X;
                               z, N_subcycle, planck_cmb, N_planck)
end


"""
    ic_int_emissivity_cmb(f_p, q, cut, bounds, V; E_X_min, E_X_max, z,
                      N_integration_steps, N_subcycle, planck_cmb, N_planck)

Integrated IC X-ray emissivity [erg s⁻¹ cm⁻³] over the photon-energy band
`[E_X_min, E_X_max]` [keV].

# Arguments
- `f_p`    : spectral norms at bin starts
- `q`      : power-law slopes per bin
- `cut`    : cutoff momentum [m_e c]
- `bounds` : N+1 bin boundaries [m_e c]

# Keyword Arguments
- `E_X_min::Real=0.1`  : lower photon-energy integration bound **[keV]**
- `E_X_max::Real=100.0`: upper photon-energy integration bound **[keV]**
                         (converted internally to erg via `keV2erg`)
- `z::Real=0.0`                    : redshift
- `N_integration_steps::Int=100`  : log-spaced photon-energy grid points over
                                    `[E_X_min, E_X_max]` (trapezoidal integration)
- `N_subcycle::Int=10`            : per-bin momentum subcycles (see `ic_source_cmb`)
- `planck_cmb::Bool=false`        : full Planck CMB spectrum vs. monochromatic ε₀
- `N_planck::Int=N_planck_default`: inner ε-integration points when `planck_cmb=true`
"""
function ic_int_emissivity_cmb(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real};
                            E_X_min::Real=0.1,
                            E_X_max::Real=100.0,
                            z::Real=0.0,
                            N_integration_steps::Int=100,
                            N_subcycle::Int=10,
                            planck_cmb::Bool=false,
                            N_planck::Int=N_planck_default)

    # E_X band supplied in keV; the emissivity machinery works in erg
    E_X_min_erg = E_X_min * keV2erg
    E_X_max_erg = E_X_max * keV2erg

    x = 10.0 .^ LinRange(log10(E_X_min_erg), log10(E_X_max_erg), N_integration_steps)
    y = [ic_emissivity_cmb(f_p, q, cut, bounds, E;
                           z, N_subcycle, planck_cmb, N_planck) for E ∈ x]

    return integrate(x, y, Trapezoidal())
end


"""
    ic_luminosity_cmb(f_p, q, cut, bounds, V; E_X_min, E_X_max, z,
                      N_integration_steps, N_subcycle, planck_cmb, N_planck)

Total IC X-ray luminosity [erg s⁻¹] for volume `V` [cm³], integrated over the
photon-energy band `[E_X_min, E_X_max]`.

# Arguments
- `f_p`    : spectral norms at bin starts
- `q`      : power-law slopes per bin
- `cut`    : cutoff momentum [m_e c]
- `bounds` : N+1 bin boundaries [m_e c]
- `V`      : emitting volume [cm³]

# Keyword Arguments
- `E_X_min::Real=0.1`  : lower photon-energy integration bound **[keV]**
- `E_X_max::Real=100.0`: upper photon-energy integration bound **[keV]**
                         (converted internally to erg via `keV2erg`)
- `z::Real=0.0`                    : redshift
- `N_integration_steps::Int=100`  : log-spaced photon-energy grid points over
                                    `[E_X_min, E_X_max]` (trapezoidal integration)
- `N_subcycle::Int=10`            : per-bin momentum subcycles (see `ic_source_cmb`)
- `planck_cmb::Bool=false`        : full Planck CMB spectrum vs. monochromatic ε₀
- `N_planck::Int=N_planck_default`: inner ε-integration points when `planck_cmb=true`
"""
function ic_luminosity_cmb(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            V::Real;
                            E_X_min::Real=0.1,
                            E_X_max::Real=100.0,
                            z::Real=0.0,
                            N_integration_steps::Int=100,
                            N_subcycle::Int=10,
                            planck_cmb::Bool=false,
                            N_planck::Int=N_planck_default)

    # integrated emissivity [erg s⁻¹ cm⁻³] over the band [E_X_min, E_X_max]
    q_ic = ic_int_emissivity_cmb(f_p, q, cut, bounds;
                            E_X_min, E_X_max, z,
                            N_integration_steps,
                            N_subcycle, planck_cmb, N_planck)

    return q_ic * V
end


"""
    ic_flux_cmb(f_p, q, cut, bounds, V, d; E_X_min, E_X_max, z,
                N_integration_steps, N_subcycle, planck_cmb, N_planck)

IC X-ray photon flux [photons s⁻¹ cm⁻²] at distance `d` [cm], integrated over the
photon-energy band `[E_X_min, E_X_max]`.

# Arguments
- `f_p`    : spectral norms at bin starts
- `q`      : power-law slopes per bin
- `cut`    : cutoff momentum [m_e c]
- `bounds` : N+1 bin boundaries [m_e c]
- `V`      : emitting volume [cm³]
- `d`      : distance to source [cm]

# Keyword Arguments
- `E_X_min::Real=0.1`  : lower photon-energy integration bound **[keV]**
- `E_X_max::Real=100.0`: upper photon-energy integration bound **[keV]**
                         (converted internally to erg via `keV2erg`)
- `z::Real=0.0`                    : redshift
- `N_integration_steps::Int=100`  : log-spaced photon-energy grid points over
                                    `[E_X_min, E_X_max]` (trapezoidal integration)
- `N_subcycle::Int=10`            : per-bin momentum subcycles (see `ic_source_cmb`)
- `planck_cmb::Bool=false`        : full Planck CMB spectrum vs. monochromatic ε₀
- `N_planck::Int=N_planck_default`: inner ε-integration points when `planck_cmb=true`
"""
function ic_flux_cmb(f_p::Vector{<:Real},
                     q::Vector{<:Real},
                     cut::Real,
                     bounds::Vector{<:Real},
                     V::Real,
                     d::Real;
                     E_X_min::Real=0.1,
                     E_X_max::Real=100.0,
                     z::Real=0.0,
                     N_integration_steps::Int=100,
                     N_subcycle::Int=10,
                     planck_cmb::Bool=false,
                     N_planck::Int=N_planck_default)

    # integrated emissivity [erg s⁻¹ cm⁻³] over the band [E_X_min, E_X_max]
    q_ic = ic_int_emissivity_cmb(f_p, q, cut, bounds;
                            E_X_min, E_X_max, z,
                            N_integration_steps,
                            N_subcycle, planck_cmb, N_planck)

    return q_ic * V / (4π * d^2)
end
