"""
    γ_emissivity_per_bin(f_p_start::Real, p_start::Real, 
                         f_p_mid::Real, p_mid::Real, 
                         f_p_end::Real, p_end::Real,
                         Eγ::Real)

Calculate the emissivity contained in a spectral bin defined by its start, mid and end values for a given photon energy `Eγ`.
"""
function γ_source_per_bin_K14(f_p_start::Real, p_start::Real, 
                          f_p_mid::Real, p_mid::Real, 
                          f_p_end::Real, p_end::Real,
                          Eγ::Real)


    # energy density at momentum p_start * integrated synchrotron kernel
    F_start = 4π * p_start^2 * f_p_start * dσγ_dEγ_K14(p_start, Eγ)

    # middle of bin
    F_mid = 4π * p_mid^2 * f_p_mid * dσγ_dEγ_K14(p_mid, Eγ)

    # end of bin
    F_end = 4π * p_end^2 * f_p_end * dσγ_dEγ_K14(p_end, Eγ)

    # bin width
    dp = p_end - p_start

    # store total synchrotron emissivity
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
    F_start = 4π * p_start^2 * f_p_start * dσγ_dEγ_Y18(p_start, Eγ)

    # middle of bin
    F_mid = 4π * p_mid^2 * f_p_mid * dσγ_dEγ_Y18(p_mid, Eγ)

    # end of bin
    F_end = 4π * p_end^2 * f_p_end * dσγ_dEγ_Y18(p_end, Eγ)

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

Source function of gamma-ray photons at energy `Eγ` in units of `N_photons erg^-1 s^-1 cm^-3` as given in Werhahn+21, Eq. A8.
"""
function gamma_source_pions(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            bounds::Vector{<:Real},
                            nH::Real; 
                            Eγ=1.0, xHe=0.76,
                            heavy_nuclei::Bool=false,
                            reduce_spectrum::Bool=true)

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
    γ_prefac = cL * a_nucl * nH

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
            p_start = p_min_γ
            f_p_start *= (p_min_γ / bounds[i])^(-q[i])
        end

        # construct mid in log-space
        p_mid = find_log_mid(p_start, p_end)

        # spectrum integration points
        f_p_mid = f_p_start * (p_mid / p_start)^(-q[i])
        f_p_end = f_p_start * (p_end / p_start)^(-q[i])

        if isnan(f_p_start)
            q_γ[i] = 0.0
            continue
        end

        # if T_p(p_end) < 10
        #     # calculate the emissivity of the bin
        #     q_γ[i] = γ_source_per_bin_Y18(f_p_start, p_start,
        #                                   f_p_mid, p_mid,
        #                                   f_p_end, p_end,
        #                                   Eγ )
        # else
            # calculate the emissivity of the bin
            q_γ[i] = γ_source_per_bin_K14(f_p_start, p_start,
                                    f_p_mid, p_mid,
                                    f_p_end, p_end,
                                    Eγ )
        #end

        if !reduce_spectrum
            bin_centers[i] = p_mid
        end
    end # loop

    if reduce_spectrum
        return γ_prefac * sum(q_γ)
    else
        return bin_centers, γ_prefac .* q_γ
    end
end
