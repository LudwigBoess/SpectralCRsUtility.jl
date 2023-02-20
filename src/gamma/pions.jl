using NumericalIntegration


"""
    γ_emissivity_per_bin(f_p_start::Real, p_start::Real, 
                         f_p_mid::Real, p_mid::Real, 
                         f_p_end::Real, p_end::Real,
                         Eγ::Real)

Calculate the emissivity contained in a spectral bin defined by its start, mid and end values for a given photon energy `Eγ`.
"""
function γ_source_per_bin_K14(f_p_start::Real, p_start::Real, 
                          p_end::Real, q::Real,
                          Eγ::Real)


    # energy density at momentum p_start * integrated synchrotron kernel
    F_start = p_start^2 * f_p_start * dσγ_dEγ_K14(T_p(p_start), Eγ)

    # middle of bin
    # construct mid in log-space
    p_mid = find_log_mid(p_start, p_end)
    # interpolate spectrum at middle of bin 
    f_p_mid = interpolate_spectrum(p_mid, f_p_start, p_start, q)
    # solve integrand
    F_mid = p_mid^2 * f_p_mid * dσγ_dEγ_K14(T_p(p_mid), Eγ)

    # end of bin
    # interpolate spectrum at end of bin 
    f_p_end = interpolate_spectrum(p_end, f_p_start, p_start, q)
    # solve integrand
    F_end = p_end^2 * f_p_end * dσγ_dEγ_K14(T_p(p_end), Eγ)

    # bin width
    dp = p_end - p_start

    # store total gamma source function
    # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
    return dp / 6 * (F_start + F_end + 4F_mid)

end

"""
    γ_emissivity_per_bin(f_p_start::Real, p_start::Real, 
                         f_p_mid::Real, p_mid::Real, 
                         f_p_end::Real, p_end::Real,
                         Eγ::Real)

Calculate the emissivity contained in a spectral bin defined by its start, mid and end values for a given photon energy `Eγ`.
"""
function γ_source_per_bin_Y18(f_p_start::Real, p_start::Real,
                            f_p_mid::Real, p_mid::Real,
                            f_p_end::Real, p_end::Real,
                            Eγ::Real)


    # energy density at momentum p_start * integrated synchrotron kernel
    F_start = p_start^2 * f_p_start * dσγ_dEγ_Y18(p_start, Eγ)

    # middle of bin
    F_mid = p_mid^2 * f_p_mid * dσγ_dEγ_Y18(p_mid, Eγ)

    # end of bin
    F_end = p_end^2 * f_p_end * dσγ_dEγ_Y18(p_end, Eγ)

    # bin width
    dp = p_end - p_start

    # store total synchrotron emissivity
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
                        reduce_spectrum::Bool=true)

Source function of gamma-ray photons at energy `Eγ` in units of `N_photons GeV^-1 s^-1 cm^-3` as given in Werhahn+21, Eq. A8.
"""
function gamma_source_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            p_min::Real, p_max::Real,
                            nH::Real, Eγ::Real=1.0;
                            Nsubcycle::Int=1,
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            reduce_spectrum::Bool=true)

    # if all norms are 0 -> q_γ = 0!
    if iszero(sum(f_p))
        return 0.0
    end

    # store number of bins 
    Nbins = length(f_p)

    # get bin size
    dp = log10(p_max / p_min) / Nbins / Nsubcycle

    # additional factor to account for interaction with heavier nuclei
    # see discussion in Werhahn+21, last part of A1.
    if heavy_nuclei
        a_nucl = 1.8 * xHe
    else
        a_nucl = 1
    end

    # prefactor to Werhahn+21, Eq. A6
    γ_prefac = 4π * m_p * cL^2 * a_nucl * nH

    # integral over pion spectrum 
    # storage array for γ emissivity
    q_γ = zeros(Nbins)

    if !reduce_spectrum
        bin_centers = Vector{Float64}(undef, Nbins)
    end

    @inbounds for i = 1:Nbins
        
        # bin integration points
        p_start_bin = p_min * 10.0^((i - 1) * dp * Nsubcycle)

        if p_start_bin > cut
            q_γ[i:end] .= 0.0
            break
        end

        p_end = p_min * 10.0^(i * dp * Nsubcycle)
        # check if bin is only partially filled
        if p_end > cut
            p_end = cut
        end

        if !reduce_spectrum
            bin_centers[i] = find_log_mid(p_start_bin, p_end)
        end
        
        # if the end of the bin does not contribute to the
        # emission we can skip the bin!
        if p_end < p_min_γ
            q_γ[i] = 0.0
            continue
        end

        if isnan(f_p[i])
            q_γ[i] = 0.0
            continue
        end

        # sub-cycle integration for increased accuracy
        for j = 1:Nsubcycle

            # start of (sub) bin
            p_start = p_min * 10.0^(((i - 1) * Nsubcycle + j - 1) * dp)

            # end of (sub) bin
            p_end   = p_min * 10.0^((i * Nsubcycle + j) * dp)
        
            # if the end of the bin does not contribute to the
            # emission we can skip the bin!
            if p_end < p_min_γ
                q_γ[i] = 0.0
                continue
            end

            # norm is defined at start of (sub) bin
            f_p_start = interpolate_spectrum(p_start, f_p[i], p_start_bin, q[i])

            # if p_min_γ is in the center of the bin, interpolate norm
            if p_start < p_min_γ
                # interolate spectrum to new start point
                f_p_start = interpolate_spectrum(p_min_γ, f_p[i], p_start_bin, q[i])
                # reset start point
                p_start = p_min_γ
            end

            
            # calculate the emissivity of the bin
            q_γ[i] += γ_source_per_bin_K14(f_p_start, p_start,
                                    p_end, q[i],
                                    Eγ )
        
        end # loop j
    end # loop i 

    if reduce_spectrum
        return γ_prefac * sum(q_γ)
    else
        return bin_centers, γ_prefac .* q_γ
    end
end


"""
    gamma_source_pions( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        bounds::Vector{<:Real},
                        nH::Real; 
                        Eγ=1.0, xHe=0.76,
                        reduce_spectrum::Bool=true)

Emissivity of gamma-ray photons at energy `Eγ` in units of `N_photons s^-1 cm^-3` as given in Werhahn+21, Eq. A2.
"""
function gamma_emissivity_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real, Eγ::Real=1.0;
                            xHe=0.76,
                            heavy_nuclei::Bool=false)

    return Eγ * gamma_source_pions(f_p, q, cut, bounds, nH, Eγ; xHe, heavy_nuclei )
end


"""
    gamma_luminosity_pions( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        bounds::Vector{<:Real},
                        nH::Real, V::Real;
                        Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                        xHe=0.76,
                        heavy_nuclei::Bool=false )

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
                            N_integration_steps::Int=100)

    x = 10.0 .^ LinRange(log10(Eγ_min), log10(Eγ_max), N_integration_steps)
    y = [gamma_emissivity_pions(f_p, q, cut, bounds, nH, Eγ; xHe, heavy_nuclei) for Eγ ∈ x]

    integral = integrate(x, y, Trapezoidal())

    return integral * V
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
                     N_integration_steps::Int=100)

Emissivity of gamma-ray photons at energy `Eγ` in units of `N_photons s^-1 cm^-3` as given in Werhahn+21, Eq. A2.
"""
function gamma_flux_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real, V::Real, d::Real;
                            Eγ_min::Real=0.2, Eγ_max::Real=300.0,
                            xHe=0.76,
                            heavy_nuclei::Bool=false,
                            N_integration_steps::Int=100)

    x = 10.0 .^ LinRange(log10(Eγ_min), log10(Eγ_max), N_integration_steps)
    y = [gamma_source_pions(f_p, q, cut, bounds, nH, Eγ; xHe, heavy_nuclei) for Eγ ∈ x]

    integral = integrate(x, y, Trapezoidal())

    return integral * V / (4π * d^2)
end