"""
    Parametrisation by Kafexhiu et. al. (2014)
    https://ui.adsabs.harvard.edu/abs/2014PhRvD..90l3014K/abstract
"""


"""
    𝒴_γ(E_γ)

Helper function as in Werhahn+21 Eq. A14
"""
𝒴_γ(Eγ) = Eγ + E_π0^2/4Eγ

"""
    𝒳_γ(E_γ)

Helper function as in Werhahn+21, Eq. A16.
"""
𝒳_γ(Tp, Eγ) = (𝒴_γ(Eγ) - E_π0) / (𝒴_γ(E_γ_max(Tp)) - E_π0)


"""
    α_K14(Tp)

Fit value from Kafexhiu+14, Tab. 5
"""
function α_K14(Tp)
    if Tp < 20
        return 1
    else
        return 1//2
    end
end


"""
    κ_K14(Tp)

Fit function for Kafexhiu+14, Eq. 14.
"""
κ_K14(Tp) = 3.29 - 0.2 * √((E_p0 / Tp)^3) 


"""
    μ_K14(Tp)

Fit function for Kafexhiu+14, Eq. 11 parameters.
"""
function μ_K14(Tp)
    q = (Tp - 1) / E_p0
    sq = √(√(q))
    return 5 / 4 * sq^2 * sq^3 * exp(-5 / 4 * q)
end



"""
    β_K14(Tp)

Fit value from Kafexhiu+14, Tab. 5
"""
function β_K14(Tp)
    if Tp < 1
        return κ_K14(Tp)
    elseif 1 < Tp ≤ 4
        return μ_K14(Tp) + 2.45
    elseif 4 < Tp ≤ 20
        return 3/2 * μ_K14(Tp) + 4.95
    elseif 20 < Tp ≤ 100
        return 4.2
    else
        return 4.9
    end
end


"""
    γ_K14(Tp)

Fit value from Kafexhiu+14, Tab. 5
"""
function γ_K14(Tp)
    if Tp < 1
        return 0
    elseif 1 < Tp ≤ 4
        return μ_K14(Tp) + 1.45
    elseif 4 < Tp ≤ 20
        return μ_K14(Tp) + 1.5
    else
        return 1
    end
end


"""
    𝒞(p)

Fit parameter from Kafexhiu+14, Eq. 11 using `Geant4` data.
"""
𝒞(Tp) = 3E_π0 / 𝒴_γ(E_γ_max(Tp))


"""
    F_K14(Tp, Eγ)

Fit function to represent cross section (Kafexhiu+14, Eq. 11)
"""
function F_K14(Tp, Eγ)

    X = 𝒳_γ(Tp, Eγ)

    if 0 ≤ X < 1
        return (1 - X^(α_K14(Tp)))^(β_K14(Tp)) / (1 + X / 𝒞(Tp))^(γ_K14(Tp))
    else
        return 0
    end
end


"""
    A_max(p)

Peak value for gamma ray spectrum following Kafexhiu+14, Eq. 12.
"""
function A_max(Tp)

    if Tp < Tp_th
        return 0
    end

    # proton kinetic energy in [GeV]
    if Tp_th ≤ Tp < 1
        A = 5.9 * σ_π0_K14(Tp) / E_π_max_LAB(Tp)
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

        θp = Tp / E_p0
        A = b1 / θp^b2 * exp(b3 * log(θp)^2) * σ_π0_K14(Tp) / E_p0
    end

    return A
end


"""
    dσγ_dEγ_K14(Tp, Eγ)

Differential gamma ray cross-section for pion-decay as a function of proton kinetic energy `Tp` and energy of resulting γ-ray `Eγ`.
See Kafexhiu+14, Eq. 8 in units of [cm^2/GeV].
"""
dσγ_dEγ_K14(Tp, Eγ) = A_max(Tp) * F_K14(Tp, Eγ)