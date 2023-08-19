
global const slope_soft   = 1.e-6
global const cL           = 2.9979e10
global const qe           = 4.8032e-10
global const m_p          = 1.6726e-24
global const m_e          = 9.10953e-28
const global erg2eV = 6.242e+11

global const γ            = 5.0/3.0

# synchrotron emission
global const mJy_factor   = 1.e26
# prefactor to the emissivity from Longair (2011), Eq. 8.128
# -> includes 4π factor from f(p) integral and (m_e * cL)^-1 from CR norm units
const global j_ν_prefac_c = 4π * √(3) * qe^3 / (m_e^2 * cL^3)

# Critican frequency from Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 20 
#  -> converted to dimensionless momentum
const global C_crit_p = 3qe / (4π * m_e * cL) 

# gamma emission
global const E_p0           = 938.272088e-3             # proton mass in [GeV]
global const E_π0           = 134.9769e-3               # π^0 rest-mass in [GeV]
global const Tp_th          = 2E_π0 + E_π0^2 / (2E_p0)  # proton threshold energy for π^0 production [GeV]
global const p_min_γ        = 1.22 / E_p0               # minimum momentum for pion production
global const GeVtoerg       = 1.60217733e-3
# one c_light cancels out!
const global γ_prefac_const = 4π / m_p #* GeVtoerg