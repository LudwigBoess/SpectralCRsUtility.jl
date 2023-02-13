"""
    T_p(p)

Kinetic proton energy in [GeV].
"""
T_p(p) = p * E_p0
#T_p(p) = √(p^2 / E_p0^2 + E_p0^2) - E_p0

"""
    θ_p(Tp)

Ratio between proton kinetic energy and rest energy.
"""
θ_p(Tp) = Tp / E_p0


"""
    𝓈(p)

Squared center of mass energy as a function of dimensionless proton momentum (See Werhahn+21, after Eq. A11).
"""
𝓈(p) = 2E_p0 * (T_p(p) + 2E_p0)


"""
    E_π_CM(p)

Total energy of the pion in the centre-of-mass system (see Werhahn+21, Eq. A10) as a function of dimensionless proton momentum.
"""
E_π_CM(p) = (𝓈(p) - (2E_p0)^2 - E_π0^2) / (2 * √(𝓈(p)))

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
E_π_max_LAB(p) = γ_CM(p) * (E_π_CM(p) + cL * P_π_CM(p) * β_CM(p))


"""
    β_CM(p)

Velocity in the center of mass system (see Werhahn+21, after Eq. A12) as a function of dimensionless proton momentum.
"""
β_CM(p) = √(1 - 1 / γ_CM(p)^2)


"""
    β_π_LAB(p)

Pion maximum velocity in the lab frame (see Kafexhio+14, after Eq. 10)
"""
β_π_LAB(p) = √(1 - (E_π0 / E_π_max_LAB(p))^2)

"""
    γ_CM(p)

Lorentz factor in the center of mass system (see Werhahn+21, Eq. A12) as a function of dimensionless proton momentum.
"""
γ_CM(p) = (T_p(p) + 2E_p0) / √(𝓈(p))


"""
    γ_π_LAB(p)

Lorentz factor of the pion in the lab system (see Kafexhiu+14, Eq. 10) as a function of dimensionless proton momentum.
"""
γ_π_LAB(p) = E_π_max_LAB(p) / E_π0


"""
    E_γ_min(p)

Minimum γ-ray energy from a pion at momentum `p`.
"""
E_γ_min(p) = E_π0 / 2 * γ_π_LAB(p) * (1 - β_π_LAB(p))

"""
    E_γ_max(p)

Maximum γ-ray energy from a pion at momentum `p`.
"""
E_γ_max(p) = E_π0 / 2 * γ_π_LAB(p) * (1 + β_π_LAB(p))

