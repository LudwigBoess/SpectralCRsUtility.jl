"""
    Werhahn+21 taken from Appendix A1
"""

"""
    f_γ_π(E_π)

Green's function for neutral pion decay (see Werhan+21, Eq. A8)
"""
f_γ_π(E_π) = 1 / √(E_π^2 - E_π0^2)


"""
    T_p(p)

Kinetic proton energy in [GeV].
"""
T_p(p) = p * E_p0

"""
    s(p)

Squared center of mass energy as a function of dimensionless proton momentum (See Werhahn+21, after Eq. A11).
"""
s(p) = 2E_p0 * (T_p(p) + 2E_p0)

"""
    γ_CM(p)

Lorentz factor in the center of mass system (see Werhahn+21, Eq. A12) as a function of dimensionless proton momentum.
"""
γ_CM(p) = (T_p(p) + 2E_p0) / √(s(p))

"""
    β_CM(p)

Velocity in the center of mass system (see Werhahn+21, after Eq. A12) as a function of dimensionless proton momentum.
"""
β_CM(p) = √(1 - 1 / γ_CM(p)^2)

"""
    E_π_CM(p)

Total energy of the pion in the centre-of-mass system (see Werhahn+21, Eq. A10) as a function of dimensionless proton momentum.
"""
E_π_CM(p) = (s(p) - (2E_p0)^2 - E_π0^2) / (2*√(s(p)))

"""
    P_π_CM(p)

Total momentum of the pion in the centre-of-mass system (see Werhahn+21, Eq. A11) as a function of dimensionless proton momentum.
"""
P_π_CM(p) = √(E_π_CM(p)^2 - E_π0^2) / cL

"""
    E_π_min(E)

Minimum allowed energy for the created pion given in the lab frame (see Werhahn+21, Eq. A17) as a function of photon energy.
"""
E_π_min(E) = max(E_π0, E + E_π0^2 / 4E)

"""
    E_π_max_LAB(p)

Maximum allowed energy for the created pion given in the lab frame (see Werhahn+21, Eq. A13) as a function of dimensionless proton momentum.
"""
E_π_max_LAB(p) = γ_CM(p) * ( E_π_CM(p) + cL * P_π_CM(p) * β_CM(p))

"""
    β_π_LAB(p)

Pion maximum velocity in the lab frame (see Kafexhio+14, after Eq. 10)
"""
β_π_LAB(p) = √(1 - (E_π0/E_π_max_LAB(p))^2)


"""
    E_γ_min(p)

Minimum γ-ray energy from a pion at momentum `p`.
"""
E_γ_min(p) = E_π0 / 2 * E_π_max_LAB(p) * (1 - β_π_LAB(p))

"""
    E_γ_max(p)

Maximum γ-ray energy from a pion at momentum `p`.
"""
E_γ_max(p) = E_π0 / 2 * E_π_max_LAB(p) * (1 + β_π_LAB(p))   


E_γ_max(1.0)


"""
    𝓃_π(p) 

Pion average yield (without ϵ!) as in Yang+18, after Eq. 5
"""
function 𝓃_π(p) 
    𝓌 = √(s(p)) / m_p
    F = (𝓌 - 2)^(3/4) / 𝓌^(1/4)
    return 0.78F - 0.5 
end

"""
    𝓃_π0(p)

Pion average yield for ``π^0`` as in Yang+18, after Eq. 5
"""
𝓃_π0(p) = 𝓃_π(p) + 1/3

"""
    𝓃𝓃_π_minus_π0(p) 

Pion average yield for ``π^-`` as in Yang+18, after Eq. 5
"""
𝓃_π_minus(p) = 𝓃_π(p)

"""
    𝓃_π_plus(p) 

Pion average yield for ``π^+`` as in Yang+18, after Eq. 5
"""
𝓃_π_plus(p)  = 𝓃_π(p) + 2 / 3


"""
    σ_pp_inel(Tp)

Inelastic cross-section for pion production as a function of proton kinetic energy (see Kafexhio+14, Eq. 1).
Caution: Werhahn+21, Eq. A20 has a typo!
"""
σ_pp_inel(Tp) = ( 30.7 - 0.96*log(Tp/Tp_th) + 0.18*log(log(Tp/Tp_th)) ) * 
                (1 - (Tp_th/Tp)^1.9)^3 * 1.e-27


"""
    σ_π0(p)

Cross-section for neutral pion decay.
"""
σ_π0(p) = σ_pp_inel(T_p(p)) * 𝓃_π0(p)

"""
    𝒴(E_γ)

Helper function as in Werhahn+21 Eq. A14
"""
𝒴_γ(E_γ) = E_γ + E_π0/4E_γ

"""
    𝒳_γ(E_γ)

Helper function as in Werhahn+21, Eq. A16.
"""
𝒳_γ(E_γ, p) = (𝒴_γ(E_γ) - E_π0) / (𝒴_γ(E_γ_max(p)) - E_π0)


p2 = 20.0
𝒳_γ(E_γ_max(p2), 1.1p2)


"""
    α(Tp)

Fit value from Kafexhio+14, Tab. 5
"""
function α_fit(Tp)
    if Tp < 20
        return 1
    else
        return 1/2
    end
end

"""
    κ(Tp)

Fit function for Kafexhio+14, Eq. 11 parameters.
"""
κ_fit(Tp) = 3.29 - 1 / 5 * √(E_p0 / Tp)^3

function μ_fit(Tp)
    q = (Tp - 1) / E_p0
    return 5 / 4 * √(√(q))^5 * exp(-5 / 4 * q)
end

"""
    β(Tp)

Fit value from Kafexhio+14, Tab. 5
"""
function β_fit(Tp)
    if Tp < 1
        return κ_fit(Tp)
    elseif 1 < Tp ≤ 4
        return μ_fit(Tp) + 2.45
    elseif 4 < Tp ≤ 20
        return 3/2 * μ_fit(Tp) + 4.95
    elseif 20 < Tp ≤ 100
        return 4.2
    else
        return 4.9
    end
end

"""
    γ_fit(Tp)

Fit value from Kafexhio+14, Tab. 5
"""
function γ_fit(Tp)
    if Tp < 1
        return 0
    elseif 1 < Tp ≤ 4
        return μ_fit(Tp) + 1.45
    elseif 4 < Tp ≤ 20
        return μ_fit(Tp) + 1.5
    else
        return 1
    end
end


"""
    𝒞(p)

Fit parameter from Kafexhio+14, Eq. 11 using `Geant4` data.
"""
𝒞(p) = 3E_π0 * 𝒴_γ(E_γ_max(p))




"""
    F(Tp, Eγ)

Fit function to represent cross section (Kafexhio+14, Eq. 11)
"""
function F(p, E_γ) 

    Tp = T_p(p)

    X = 𝒳_γ(E_γ, p)

    if 0 ≤ X < 1
        return (1 - X^(α_fit(Tp)))^(β_fit(Tp)) / (1 + X / 𝒞(p))^(γ_fit(Tp))
    else
        return 0
    end
end


"""
    A_max(p)

Peak value for gamma ray spectrum following Kafexhio+14, Eq. 12.
"""
function A_max(p)

    Tp = T_p(p)

    # proton kinetic energy in [GeV]
    if Tp_th ≤ Tp < 1
        A = 5.9 * σ_π0(p) / E_π_max_LAB(p)
    else
        # we use Geant4 results from Tab. VII.
        if 1 ≤ Tp < 5
            b1 = 9.53 
            b2 = 0.52 
            b3 = 0.054
        else
            b1 = 9.13 
            b2 = 0.35 
            b3 = 9.7e-3
        end
        θ_p = Tp / E_p0
        A = b1 * θ_p^(-b2) * exp(b3 * log(log(θ_p)) ) * σ_π0(p) / E_p0
    end

    return A
end


"""
    dσ_dEγ(p, Eγ)

Differential cross-section for pion-decay as a function of dimensionless proton momentum `p` and energy of resulting γ-ray `Eγ`.
See Kafexhio+14, Eq. 8
"""
dσπ_dEγ(p, Eγ) = A_max(p) * F(T_p(p), Eγ)




"""
    dσγ_dEγ_trapez(p::Real, Eγ::Real, E_steps::Integer=50)

Pitch angle integration in Donnert+16, Eq. 17.
"""
function dσγ_dEγ_trapez(p::Real, Eγ::Real, E_steps::Integer=50)

    Eπ_min = E_π_min(Eγ)
    Eπ_max = E_π_max_LAB(p)

    dEπ = (Eπ_max - Eπ_min) / E_steps
    Eπ = Eπ_min

    # first half step: sin(0) = 0
    dσγ_dEγs = 0.5 * dσπ_dEγ(p, Eγ) * f_γ_π(Eπ) 

    # actual integration
    @inbounds for i ∈ 1:E_steps-1
        Eπ += dEπ
        dσγ_dEγs += dσπ_dEγ(p, Eγ) * f_γ_π(Eπ)
    end

    # last step
    Eπ += dEπ
    dσγ_dEγs += 0.5 * dσπ_dEγ(p, Eγ) * f_γ_π(Eπ)

    # multiply by step length
    dσγ_dEγs *= dEπ

    # factor 2 comes from pions decaying into 2 photons
    return 2dσγ_dEγs
end


"""
    γ_emissivity_per_bin(f_p_start::Real, p_start::Real, 
                         f_p_mid::Real, p_mid::Real, 
                         f_p_end::Real, p_end::Real,
                         Eγ::Real)

Calculate the emissivity contained in a spectral bin defined by its start, mid and end values for a given photon energy `Eγ`.
"""
function γ_source_per_bin(f_p_start::Real, p_start::Real, 
                                f_p_mid::Real, p_mid::Real, 
                                f_p_end::Real, p_end::Real,
                                Eγ::Real)


    # energy density at momentum p_start * integrated synchrotron kernel
    F_start = 4π * p_start^2 * f_p_start * dσγ_dEγ_trapez(p_start, Eγ)

    # middle of bin
    F_mid = 4π * p_mid^2 * f_p_mid * dσγ_dEγ_trapez(p_mid, Eγ)

    # end of bin
    F_end = 4π * p_end^2 * f_p_end * dσγ_dEγ_trapez(p_end, Eγ)

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

        # calculate the emissivity of the bin
        q_γ[i] = γ_source_per_bin(f_p_start, p_start,
                                  f_p_mid, p_mid,
                                  f_p_end, p_end,
                                  Eγ )

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
