# Thomson cross-section [cm²]
const σ_T = 6.6524e-25

# Boltzmann constant [erg/K]
const k_B = 1.3807e-16

# Stefan-Boltzmann constant [erg s⁻¹ cm⁻² K⁻⁴]
const σ_SB = 5.6704e-5

# CMB temperature today [K]
const T_CMB0 = 2.7255

# Mean CMB photon energy factor: ε₀ = ic_mean_photon_factor × k_B × T
const ic_mean_photon_factor = 2.7

# Prefactor: 3π σ_T / m_e
# Absorbs 4π from f(p) solid-angle integral, 3/4 from Thomson IC formula,
# and 1/(m_e cL) from CR norm units (same convention as j_ν_prefac_c).
const ic_prefac_c = 3π * σ_T / m_e

# Planck spectrum prefactor: 8π / (h c)³
# n_CMB(ε, T) = ic_planck_prefac · ε² / (exp(ε/k_B T) - 1)  [photons cm⁻³ erg⁻¹]
const ic_planck_prefac = 8π / (h_planck * cL)^3

# Default number of points for the inner ε-integration of the Planck spectrum
const N_planck_default = 64
