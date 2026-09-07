"""
    T_p(p)

Kinetic proton energy in [GeV].
"""
T_p(p) = E_p0 * (√(1 + p^2) - 1)


"""
    θ_p(Tp)

Ratio between proton kinetic energy and rest energy.
"""
θ_p(Tp) = Tp / E_p0


"""
    𝓈(Tp)

Squared center of mass energy as a function of proton kinetic energy (See Werhahn+21, after Eq. A11).
"""
𝓈(Tp) = 2E_p0 * (Tp + 2E_p0)


"""
    E_π_CM(Tp)

Total energy of the pion in the centre-of-mass system (see Werhahn+21, Eq. A10) as a function of proton kinetic energy.
"""
E_π_CM(Tp) = (𝓈(Tp) - (2E_p0)^2 + E_π0^2) / (2 * √(𝓈(Tp)))



"""
    P_π_CM(p)

Total momentum of the pion in the centre-of-mass system (see Werhahn+21, Eq. A11) as a function of proton kinetic energy.

The pion is produced at rest in the centre-of-mass system at the production threshold, where `p_min_γ`
starts the spectral integral. The `max` guards the exact cancellation there against round-off pushing
the radicand below zero.
"""
P_π_CM(Tp) = √(max(0.0, E_π_CM(Tp)^2 - E_π0^2))


"""
    E_π_min(E)

Minimum allowed energy for the created pion given in the lab frame (see Werhahn+21, Eq. A17) as a function of photon energy.
"""
E_π_min(E) = max(E_π0, E + E_π0^2 / 4E)


"""
    E_π_max_LAB(p)

Maximum allowed energy for the created pion given in the lab frame (see Werhahn+21, Eq. A13) as a function of proton kinetic energy.
"""
E_π_max_LAB(Tp) = γ_CM(Tp) * (E_π_CM(Tp) + P_π_CM(Tp) * β_CM(Tp))


"""
    β_CM(Tp)

Velocity in the center of mass system (see Werhahn+21, after Eq. A12) as a function of proton kinetic energy.
"""
β_CM(Tp) = √(1 - 1 / γ_CM(Tp)^2)


"""
    β_π_LAB(Tp)

Pion maximum velocity in the lab frame (see Kafexhio+14, after Eq. 10)
"""
β_π_LAB(Tp) = √(1 - (E_π0 / E_π_max_LAB(Tp))^2)


"""
    γ_CM(Tp)

Lorentz factor in the center of mass system (see Werhahn+21, Eq. A12) as a function of proton kinetic energy.
"""
γ_CM(Tp) = (Tp + 2E_p0) / √(𝓈(Tp))


"""
    γ_π_LAB(Tp)

Lorentz factor of the pion in the lab system (see Kafexhiu+14, Eq. 10) as a function of proton kinetic energy.
"""
γ_π_LAB(Tp) = E_π_max_LAB(Tp) / E_π0


"""
    E_γ_min(Tp)

Minimum γ-ray energy from a pion with kinetic energy `Tp`.
"""
E_γ_min(Tp) = E_π0 / 2 * γ_π_LAB(Tp) * (1 - β_π_LAB(Tp))


"""
    E_γ_max(Tp)

Maximum γ-ray energy from a pion with kinetic energy `Tp`.
"""
E_γ_max(Tp) = E_π0 / 2 * γ_π_LAB(Tp) * (1 + β_π_LAB(Tp))
