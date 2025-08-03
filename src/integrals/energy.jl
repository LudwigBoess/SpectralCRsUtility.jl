

@inline function energy_integral(bound_low::Real, bound_up::Real,
                                 norm::Real, slope::Real, ρ::Real)
    
    # energy integral (eq.21 M01)
    en = 4π * cL * norm * bound_low^4 /  ρ
    energy = en * ( (bound_up/bound_low)^(4 - slope ) - 1 ) / ( 4 - slope )

    if ( 4 - slope_soft ) < slope < ( 4 + slope_soft )
        slope_var = (slope - 4)/slope_soft
        energy2 = en * log10(bound_up/bound_low)
        if slope_var != 0
            energy = energy * slope_var + energy2 * ( 1 - slope_var )
        else
            energy = energy2
        end
    end

    return energy
end

energy_integrand(p::T, q::T) where T<: Real = (√(p^2 + 1) - 1) * p^(2 - q)

@inline function energy_integral_exact(bound_low::Real, bound_up::Real, slope::Real)
    
    # energy integral without relativistic limit 
    n_intervals = max(6, round(Int, abs(log10(bound_up/bound_low)) * 12))
    ds = log(bound_up/bound_low) / n_intervals # equal stepsize
    exp_ds = exp(ds)
    exp_2ds = exp_ds * exp_ds

    # furst summand 
    p = bound_low
    e_sum = p

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