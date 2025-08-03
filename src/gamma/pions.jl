using NumericalIntegration


"""
    γ_emissivity_per_bin(f_p_start::Real, p_start::Real, 
                         f_p_mid::Real, p_mid::Real, 
                         f_p_end::Real, p_end::Real,
                         Eγ::Real)

Calculate the emissivity contained in a spectral bin defined by its start, mid and end values for a given photon energy `Eγ`.
"""
function γ_source_per_bin(f_p_start::Real, p_start::Real, 
                          p_end::Real, q::Real,
                          Eγ::Real, dσγ_dEγ::Function)


    # energy density at momentum p_start * integrated synchrotron kernel
    F_start = p_start^2 * f_p_start * dσγ_dEγ(T_p(p_start), Eγ)

    # middle of bin
    # construct mid in log-space
    p_mid = find_log_mid(p_start, p_end)
    # interpolate spectrum at middle of bin 
    f_p_mid = interpolate_spectrum(p_mid, f_p_start, p_start, q)
    # solve integrand
    F_mid = p_mid^2 * f_p_mid * dσγ_dEγ(T_p(p_mid), Eγ)

    # end of bin
    # interpolate spectrum at end of bin 
    f_p_end = interpolate_spectrum(p_end, f_p_start, p_start, q)
    # solve integrand
    F_end = p_end^2 * f_p_end * dσγ_dEγ(T_p(p_end), Eγ)

    # bin width
    dp = p_end - p_start

    # store total gamma source function
    # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
    return dp / 6 * (F_start + F_end + 4F_mid)

end


"""
    gamma_source_pions( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        bounds::Vector{<:Real},
                        nH::Real; 
                        Eγ=1.0, xHe=0.76,
                        reduce_spectrum::Bool=true,
                        N_subcycle::Int=10)

Source function of gamma-ray photons at energy `Eγ` in units of `N_photons GeV^-1 s^-1 cm^-3` as given in Werhahn+21, Eq. A8.
"""
function gamma_source_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real, Eγ::Real=1.0;
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            reduce_spectrum::Bool=true,
                            N_subcycle::Int=10)

    # if all norms are 0 -> q_γ = 0!
    if iszero(sum(f_p))
        return 0.0
    end

    # store number of bins 
    Nbins = length(f_p)

    # additional factor to account for interaction with heavier nuclei
    # see discussion in Werhahn+21, last part of A1.
    if heavy_nuclei
        a_nucl = 1.8 * xHe
    else
        a_nucl = 1
    end

    # prefactor to Werhahn+21, Eq. A6
    γ_prefac = γ_prefac_const * a_nucl * nH

    # integral over pion spectrum 
    # storage array for γ emissivity
    q_γ = Vector{Float64}(undef, Nbins)

    if !reduce_spectrum
        bin_centers = Vector{Float64}(undef, Nbins)
    end

    @inbounds for i = 1:Nbins

        # bin integration points
        p_start = bounds[i]

        if p_start > cut
            q_γ[i:end] .= 0.0
            break
        end

        p_end = bounds[i+1]
        # check if bin is only partially filled
        if p_end > cut
            p_end = cut
        end

        if !reduce_spectrum
            bin_centers[i] = find_log_mid(p_start, p_end)
        end
        
        # if the end of the bin does not contribute to the
        # emission we can skip the bin!
        if p_end < p_min_γ
            q_γ[i] = 0.0
            continue
        end

        # norm is defined at start of bin
        f_p_start = f_p[i]

        # if p_min_γ is in the center of the bin, interpolate norm
        if p_start < p_min_γ
            # interolate spectrum to new start point
            f_p_start = interpolate_spectrum(p_min_γ, f_p_start, p_start, q[i])
            # reset start point
            p_start = p_min_γ
        end

        if isnan(f_p_start)
            q_γ[i] = 0.0
            continue
        end

        if N_subcycle > 1
            q_γ[i] = 0.0

            p_subcycle = 10.0.^LinRange(log10(p_start), log10(p_end), N_subcycle + 1)
            # loop over subcycled bins 
            for j = 1:N_subcycle
                f_p_subcycle = interpolate_spectrum(p_subcycle[j], f_p_start, p_start, q[i])

                q_γ[i] += γ_source_per_bin(f_p_subcycle, p_subcycle[j],
                                        p_subcycle[j+1], q[i],
                                        Eγ, dσγ_dEγ_K14)
            end
        else

        # if T_p(p_end) < 10.0
        #     # calculate the emissivity of the bin
        #     q_γ[i] = γ_source_per_bin(f_p_start, p_start,
        #                                 p_end, q[i],
        #                                 Eγ, dσγ_dEγ_Y18)
        # else
            # calculate the emissivity of the bin
            q_γ[i] = γ_source_per_bin(f_p_start, p_start,
                                      p_end, q[i],
                                      Eγ, dσγ_dEγ_K14)
        # end
        end # if subcycle
        
    end # loop

    if reduce_spectrum
        return γ_prefac * sum(q_γ)
    else
        return bin_centers, γ_prefac .* q_γ
    end
end


"""
    gamma_emissivity_pions( f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real, Eγ::Real=1.0;
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            N_subcycle::Int=10)

Emissivity of gamma-ray photons at energy `Eγ` in units of `N_photons s^-1 cm^-3` as given in Werhahn+21, Eq. A2.
"""
function gamma_emissivity_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real, Eγ::Real=1.0;
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            N_subcycle::Int=10)

    return Eγ * gamma_source_pions(f_p, q, cut, bounds, nH, Eγ; xHe, heavy_nuclei, N_subcycle )
end


"""
    gamma_luminosity_pions( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        bounds::Vector{<:Real},
                        nH::Real, V::Real;
                        Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                        xHe=0.76,
                        heavy_nuclei::Bool=false,
                        N_subcycle::Int=10)

Total gamma-ray luminosity of in units of `erg s^-1` as given in Werhahn+21, Eq. A4.
"""
function gamma_luminosity_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real, V::Real;
                            Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            N_integration_steps::Int=100,
                            N_subcycle::Int=10)

    x = 10.0 .^ LinRange(log10(Eγ_min), log10(Eγ_max), N_integration_steps)
    y = [gamma_emissivity_pions(f_p, q, cut, bounds, nH, Eγ; xHe, heavy_nuclei, N_subcycle) for Eγ ∈ x]

    integral = integrate(x, y, Trapezoidal())

    return integral * V * GeVtoerg
end

"""
    gamma_flux_pions(f_p::Vector{<:Real},
                     q::Vector{<:Real},
                     cut::Real,
                     bounds::Vector{<:Real},
                     nH::Real, V::Real, d::Real;
                     Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                     xHe=0.76,
                     heavy_nuclei::Bool=false,
                     N_integration_steps::Int=100,
                     N_subcycle::Int=10)

Flux of gamma-ray photons at energy `Eγ` in units of `N_photons s^-1 cm^-2` as given in Werhahn+21, Eq. A2.
"""
function gamma_flux_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real, V::Real, d::Real;
                            Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            N_integration_steps::Int=100,
                            N_subcycle::Int=10)

    x = 10.0 .^ LinRange(log10(Eγ_min), log10(Eγ_max), N_integration_steps)
    y = [gamma_source_pions(f_p, q, cut, bounds, nH, Eγ; xHe, heavy_nuclei, N_subcycle) for Eγ ∈ x]

    integral = integrate(x, y, Trapezoidal())

    return integral * V / (4π * d^2)
end

"""
    gamma_flux_pions(f_p::Vector{<:Real},
                     q::Vector{<:Real},
                     cut::Real,
                     par::CRMomentumDistributionConfig,
                     nH::Real, V::Real, d::Real;
                     Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                     xHe=0.76,
                     heavy_nuclei::Bool=false,
                     N_integration_steps::Int=100,
                     N_subcycle::Int=10)

Flux of gamma-ray photons at energy `Eγ` in units of `N_photons s^-1 cm^-2` as given in Werhahn+21, Eq. A2.
"""
function gamma_flux_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            par::CRMomentumDistributionConfig,
                            nH::Real, V::Real, d::Real;
                            Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            N_integration_steps::Int=100,
                            N_subcycle::Int=10)

    # construct boundaries
    bounds = momentum_bin_boundaries(par)

    return gamma_flux_pions(f_p, q, cut, bounds, nH, V, d; 
                            Eγ_min, Eγ_max, 
                            xHe, heavy_nuclei, 
                            N_integration_steps,
                            N_subcycle)
end