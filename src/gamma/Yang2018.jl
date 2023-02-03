"""
    Parametrisation by Yang et. al. (2018)
    https://ui.adsabs.harvard.edu/abs/2018A%26A...615A.108Y/abstract
"""

"""
    Fit functions in Appendix
"""

"""
    𝒮(x)

Sigmoid function
"""
𝒮(x) = 1 / (1 + exp(-x))



"""
    x0_Y18(Tp)

Fit function from Yang+18, Eq. A4
"""
function x0_Y18(Tp)
    if Tp < 1.5
        return 0.17
    elseif 1.5 ≤ Tp < 5
        return 0.1
    elseif 5 ≤ Tp < 10
        return 0.08
    end
end


"""
    α_Y18(Tp)

Fit function from Yang+18, Eq. A5
"""
function α_Y18(Tp)
    if Tp < 1
        return 0.42 - 0.1*θ_p(Tp)
    elseif 1 ≤ Tp < 1.5
        return 0.36
    elseif 1.5 ≤ Tp < 10
        θp = θ_p(Tp)
        return 0.288 + (𝒮(18 * (θp - 1.85)) + 10 * 𝒮(1.7*(θp - 4.1))) / 40
    end
end

"""
    β_Y18(Tp)

Fit function from Yang+18, Eq. A6
"""
function β_Y18(Tp)
    if Tp < 1
        θp = θ_p(Tp)
        return 101 * θp^3 - 230 * θp^2 + 170 * θp - 30
    elseif 1 ≤ Tp < 1.5
        return 12.3
    elseif 1.5 ≤ Tp < 10
        θp = θ_p(Tp)
        return 18.5 + 10 * 𝒮(13 * (θp - 3.4))
    end
end

"""
    γ_Y18(Tp)

Fit function from Yang+18, Eq. A6
"""
function γ_Y18(Tp)
    if Tp < 1
        return 2.5 - 2θ_p(Tp)
    elseif 1 ≤ Tp < 1.5
        return 0.68
    elseif 1.5 ≤ Tp < 10
        return 1.9 * √(θ_p(Tp)) - 1.1
    end
end


"""
    𝒩_low(Tp)

Fit function from Yang+18, Eq. A3
"""
function 𝒩_low(Tp)

    x0 = x0_Y18(Tp)
    β  = β_Y18(Tp)

    return β^2 / (1 - (1 - β*x0)*exp(-β*x0))
end

"""
    f_low(x, Tp)

Fit function from Yang+18, Eq. A2
"""
f_low(x, Tp) = 𝒩_low(Tp) * x * exp(-β_Y18(Tp) * x)

"""
    𝒩_high(Tp)

Fit function from Yang+18, Eq. A3
"""
function 𝒩_high(Tp)

    x0 = x0_Y18(Tp)
    γ  = γ_Y18(Tp)
    
    return γ^2 / (exp(-γ) - (1 - γ * (1 - x0)) * exp(-γ*x0))
end


"""
    f_high(x, Tp)

Fit function from Yang+18, Eq. A2
"""
f_high(x, Tp) = 𝒩_high(Tp) * (1 - x) * exp(-γ_Y18(Tp) * x)

"""
    f_Y18(x, Tp)

Fit function from Yang+18, Eq. A1.
`x = E_π / E_π_max` is the ratio between pion energy and maximum pion energy in the lab frame.
`Tp` is the proton kinetic energy.
"""
function f_Y18(x, Tp)

    x0 = x0_Y18(Tp)
    
    if x ≤ x0
        return α_Y18(Tp) * f_low(x, Tp)
    elseif x0 < x ≤ 1
        return (1 - α_Y18(Tp)) * f_high(x, Tp)
    else
        return 0
    end
end



"""
    f_γ_π(E_π)

Green'𝓈 function for neutral pion decay (see Werhan+21, Eq. A8)
"""
f_γ_π(E_π) = 1 / √(E_π^2 - E_π0^2)


using QuadGK

"""
    dσγ_dEγ_Y18(p, Eγ)
"""
function dσγ_dEγ_Y18(p, Eγ)

    Eπ_min = E_π_min(Eγ)
    Eπ_max = E_π_max_LAB(p)

    σπ0 = σ_π0_Y18(p)

    Tp = T_p(p)

    # solve integral numerically
    integral, err = quadgk(Eπ -> f_Y18(Eπ / Eπ_max, Tp), Eπ_min, Eπ_max, rtol=1e-8)

    return 2σπ0 * Eπ_max * integral
end

# E = 1.0
# p = 2.0

# @btime dσγ_dE($E, $p)


# values = [0.5, 1.2, 2.5, 7.2]
# for value ∈ values
#     println(f_Y18(0.5, value))
# end