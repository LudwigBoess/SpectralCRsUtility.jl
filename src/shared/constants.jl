
global const slope_soft   = 1.e-6
global const c_light      = 2.9979e10
global const qe           = 4.8032e-10
global const m_p          = 1.6726e-24
global const m_e          = 9.10953e-28
global const nu_c_prefac  = 3.0 * qe / (4π * m_e^3 * c_light^5)
global const C_crit       = 3qe / ( 4π * m_e * c_light ) # Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 20 
                                                         #  -> converted to dimensionless momentum
global const γ            = 5.0/3.0

# synchrotron emission
global const mJy_factor   = 1.e26
global const erg2eV       = 6.242e+11
global const j_ν_prefac_c = √(3) * qe^3 / c_light

# gamma emission
global const mπ_c2        = 134.9769e-3 # π^0 rest-mass in [GeV]