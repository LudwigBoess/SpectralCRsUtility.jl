
"""
    energy_integral(bound_low::Real, bound_up::Real,
                    norm::Real, slope::Real, ρ::Real)

CR energy per unit mass between `bound_low` and `bound_up` for the spectrum
`f(p) = norm (p/bound_low)^(-slope)`, in the ultra-relativistic limit `T(p) = p mc²`.
"""
@inline function energy_integral(bound_low::Real, bound_up::Real,
                                 norm::Real, slope::Real, ρ::Real)

    # energy integral (eq.21 M01)
    en = 4π * cL * norm * bound_low^4 /  ρ
    energy = en * ( (bound_up/bound_low)^(4 - slope ) - 1 ) / ( 4 - slope )

    if ( 4 - slope_soft ) < slope < ( 4 + slope_soft )
        slope_var = (slope - 4)/slope_soft
        energy2 = en * log(bound_up/bound_low)
        if slope_var != 0
            energy = energy * slope_var + energy2 * ( 1 - slope_var )
        else
            energy = energy2
        end
    end

    return energy
end

"""
    energy_integrand(p::Real, q::Real)

Integrand of the CR energy for the exact kinetic energy `T(p) = (√(p² + 1) - 1) mc²`
and a spectrum `f(p) ∝ p^(-q)`.
"""
energy_integrand(p::Real, q::Real) = (√(p^2 + 1) - 1) * p^(2 - q)

"""
    energy_integral_exact(bound_low::Real, bound_up::Real, slope::Real)

Integral of [`energy_integrand`](@ref) from `bound_low` to `bound_up`, without the prefactor
`4πc norm bound_low^slope / ρ`. Simpson's rule in `ln p` with 12 steps per full decade, as
`LookupCREnergy::exact_energy_integral` of OpenGadget3.
"""
@inline function energy_integral_exact(bound_low::Real, bound_up::Real, slope::Real)

    # energy integral without relativistic limit
    n_intervals = max(6, trunc(Int, abs(log10(bound_up/bound_low))) * 12)
    ds = log(bound_up/bound_low) / n_intervals # equal stepsize
    exp_ds = exp(ds)
    exp_2ds = exp_ds * exp_ds

    # first summand
    p = bound_low
    e_sum = energy_integrand(p, slope) * p

    # weights for odd indices are 4, for even indices 2

    # odd indices
    p *= exp_ds
    for i in 1:2:n_intervals-1
        e_sum += 4 * energy_integrand(p, slope) * p
        p *= exp_2ds
    end

    # even indices
    p = bound_low * exp_2ds
    for i in 2:2:n_intervals-1
        e_sum += 2 * energy_integrand(p, slope) * p
        p *= exp_2ds
    end
    # last summand
    e_sum += energy_integrand(bound_up, slope) * bound_up

    return e_sum * ds / 3
end

"""
    power_law_integral(bound_low::Real, bound_up::Real, norm::Real, slope::Real)

Integral of `norm p^(-slope)` from `bound_low` to `bound_up`.
"""
@inline function power_law_integral(bound_low::Real, bound_up::Real, norm::Real, slope::Real)

    integral = ( bound_up^(1 - slope) - bound_low^(1 - slope) ) / ( 1 - slope )

    # blend with the logarithmic limit to avoid the division by zero at slope = 1
    if abs(slope - 1) < slope_soft
        slope_var = (slope - 1) / slope_soft
        integral2 = log(bound_up / bound_low)
        if slope_var != 0
            integral = integral * slope_var + integral2 * ( 1 - slope_var )
        else
            integral = integral2
        end
    end

    return integral * norm
end

"""
    kinetic_energy_power_law(bin_low::Real, bin_up::Real)

Power law `N p^(-s)` that represents `energy_integrand(p, 0)` in the momentum bin
`[bin_low, bin_up]`: `s` connects the values at the bin boundaries in log-log space and `N`
matches the value at the logarithmic bin centre. Returns `(N, s)`.
"""
@inline function kinetic_energy_power_law(bin_low::Real, bin_up::Real)

    s = -( log10(energy_integrand(bin_up, 0.0)) - log10(energy_integrand(bin_low, 0.0)) ) /
         ( log10(bin_up) - log10(bin_low) )

    p_mid = √(bin_low * bin_up)
    N = energy_integrand(p_mid, 0.0) * p_mid^s

    return N, s
end

"""
    energy_integral_dk(bin_low::Real, bin_up::Real,
                       bound_low::Real, bound_up::Real,
                       norm::Real, slope::Real, ρ::Real)

CR energy per unit mass between `bound_low` and `bound_up` inside the momentum bin
`[bin_low, bin_up]`, for the spectrum `f(p) = norm (p/bound_low)^(-slope)` and the exact
kinetic energy `T(p) = (√(p² + 1) - 1) mc²`. The kinetic-energy integrand is replaced by its
power law over the bin ([`kinetic_energy_power_law`](@ref)), as in
`LookupCREnergy::calculate_CR_energy` of OpenGadget3 with `DK_ENERGY_INTEGRAL`.
"""
@inline function energy_integral_dk(bin_low::Real, bin_up::Real,
                                    bound_low::Real, bound_up::Real,
                                    norm::Real, slope::Real, ρ::Real)

    N, s = kinetic_energy_power_law(bin_low, bin_up)

    return 4π * cL * norm * bound_low^slope / ρ * power_law_integral(bound_low, bound_up, N, s + slope)
end




# function calculateCREnergyInCGS(CR_N, CR_S, CR_Cut, ρ;
#                                 pmin::Real=10.0, pmax::Real=1.e7,
#                                 mc::Real=0.0, SelectBin::Int=-1,
#                                 units::GadgetPhysical=GadgetPhysical(),
#                                 verbose=true)
#     # calculates engergy per bin in cgs units.

#     if mc == 0.0
#         error("Error! No mc specified! Please provide m_particle * cL to continue!")
#         return 1
#     end

#     CR_E = 0.0
#     cL = 2.9979e10
#     ρ *= units.m_unit / units.l_unit^3  # convert rho in cgs units

#     CR_N .*= 1.e20

#     Nbins = size(CR_N, 1)
#     bin_width = log10(pmax/pmin)/Nbins
#     bin = 10.0.^collect(log10(pmin):bin_width:log10(pmax)) .* mc
#     above_cut = findall( bin .> CR_Cut/mc )
#     if length(above_cut) > 0
#         bin[above_cut] .= CR_Cut/mc
#     end

#     if SelectBin == -1
#         CR_E = zeros(length(CR_N))
#         for i = 1:length(CR_N)
#             CR_E[i] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
#         end
#     elseif length(SelectBin) == 1
#         CR_E = energy_integral(bin[SelectBin], bin[SelectBin+1], CR_N[SelectBin], CR_S[SelectBin], ρ)
#     elseif length(SelectBin) > 1
#         CR_E = zeros(length(SelectBin))
#         if i ∈ SelectBin
#             CR_E[i-SelectBin[1]+1] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
#         end
#     end
#     return CR_E
# end