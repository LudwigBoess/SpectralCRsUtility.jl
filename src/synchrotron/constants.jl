"""
    Constants for synchrotron
"""

const c_light    = 2.9979e10
const m_e        = 9.10953e-28
const q_e        = 1.602176487e-20 * c_light
const C_crit_p   = 3q_e / ( 4π * m_e * c_light ) # Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 20 
                                                        #  -> converted to dimensionless momentum
const γ          = 5.0/3.0
const mJy_factor = 1.e26       # conversion factor from [erg/cm^3/Hz/s] to mJy/cm.