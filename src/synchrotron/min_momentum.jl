
"""
    smallest_synch_bright_p(ν::Real, B::Real)

Calculate the smallest synchrotron bright moment following Donnert+16, Eq. 22
We assume x = 100
"""
function smallest_synch_bright_p(ν::Real, B::Real)
    sqrt( ν / (100.0*C_crit_p * B) )
end