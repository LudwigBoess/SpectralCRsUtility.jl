"""
    𝓃_π(p) 

Pion average yield (without ϵ!) as in Yang+18, after Eq. 5
"""
function 𝓃_π(p)
    𝓌 = √(𝓈(p)) / m_p
    F = (𝓌 - 2)^(3 / 4) / 𝓌^(1 / 4)
    return 0.78F - 0.5
end

"""
    𝓃_π0(p)

Pion average yield for ``π^0`` as in Yang+18, after Eq. 5
"""
𝓃_π0(p) = 𝓃_π(p) + 1 / 3

"""
    𝓃𝓃_π_minus_π0(p) 

Pion average yield for ``π^-`` as in Yang+18, after Eq. 5
"""
𝓃_π_minus(p) = 𝓃_π(p)

"""
    𝓃_π_plus(p) 

Pion average yield for ``π^+`` as in Yang+18, after Eq. 5
"""
𝓃_π_plus(p) = 𝓃_π(p) + 2 / 3


"""
    σ_pp_inel(Tp)

Inelastic cross-section for pion production as a function of proton kinetic energy in `[cm^2]` (see Kafexhio+14, Eq. 1).
Caution: Werhahn+21, Eq. A20 has a typo!
"""
σ_pp_inel(Tp) = (30.7 - 0.96 * log(Tp / Tp_th) + 0.18 * log(Tp / Tp_th)^2) *
                (1 - (Tp_th / Tp)^1.9)^3 * 1.e-27


"""
    σ_π0_Y18(p)

Cross-section for neutral pion decay.
"""
σ_π0_Y18(p) = σ_pp_inel(T_p(p)) * 𝓃_π0(p)


"""
    Kafexhiu+14 parameters 
"""

"""
    η_K14(p)

Fit function given in in Kafexhiu+14, Eq. 3.
"""
η_K14(s) = √((s - E_π0^2 - 4 * E_p0^2)^2 - 16 * E_π0^2 * E_p0^2) / (2E_π0 * √(s) )

# resonant mass
const M_res = 1.1883 # [GeV]
# resonant length
const Γ_res = 0.2264 # [GeV]
# γ value for Breit–Wigner distribution
const γ_BW = √(M_res^2 * (M_res^2 + Γ_res^2))
# K value for Breit–Wigner distribution
const 𝒦_BW = √(8) * M_res * Γ_res * γ_BW / (π * √(M_res^2 + γ_BW))

"""
    f_BW(p)

Unit-less relativistic Breit–Wigner distribution as a function of center of mass energy `s` (see Kafexhiu+14, Eq. 4).
"""
f_BW(s) = E_p0 * 𝒦_BW / ( ( ( √(s) - E_p0 )^2 - M_res^2)^2 + M_res^2*Γ_res^2)


"""
    σ1π(Tp)

Cross-section for ``pp -> pp2π^0`` and ``pp -> {pn,D} π^+ π^0`` as given in Kafexhiu+14, Eq. 5 in `[cm^2]`.
"""
function σ1π(Tp)
    #s = 𝓈(p)
    s = 2E_p0*(Tp + 2E_p0)
    η = η_K14(s)
    7.66e-30 * η^1.95 * (1 + η + η^5) * f_BW(s)^1.86
end


"""
    σ2π(Tp)

Cross-section for ``pp -> pp2π^0`` and ``pp -> {pn,D} π^+ π^0`` as given in Kafexhiu+14, Eq. 5 in `[cm^2]`.
"""
function σ2π(Tp)
    
    if 0.56 ≤ Tp ≤ 2
        return 5.7e-27 / (1 + exp(-9.3 * (Tp - 1.4)))
    else
        return 0
    end
end

"""
    𝓃_π0_K14(Tp)

Pion average yield for ``π^0`` as in Kafexhiu+14, Eq. 6 & 7.
"""
function 𝓃_π0_K14(Tp)

    if 1 ≤ Tp < 5
        Qp = (Tp - Tp_th)/E_p0
        return -6.0e-3 + 0.237Qp - 0.023 * Qp^2
    elseif Tp ≥ 5
        ξ = (Tp - 3) / E_p0
        # we use Geant 4 values from Tab. IV
        a1 = 0.728 
        a2 = 0.596 
        a3 = 0.491 
        a4 = 0.2503 
        a5 = 0.117
        return a1 * ξ^a4 * (1 + exp(-a2 * ξ^a5)) * (1 - exp(-a3 * √(√(ξ))))
    else
        return 0
    end
end

"""
    σ_π0_K14(p)

Cross-section for neutral pion decay.
"""
function σ_π0_K14(Tp)
    if Tp_th ≤ Tp < 2
        return σ1π(Tp) + σ2π(Tp)
        # we use Tp_tran = 1.e5, the Geant 4.10.0 value
    elseif 2 ≤ Tp ≤ 1.e6
        return 𝓃_π0_K14(Tp) * σ_pp_inel(Tp)
    else 
        return 0
    end
end