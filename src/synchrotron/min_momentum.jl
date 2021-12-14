
"""
    get_p_min(ν::Real, B::Real)

Calculate the smallest synchrotron bright moment following Donnert+16, Eq. 22
"""
function smallest_synch_bright_p(ν::Real, B::Real)
    sqrt( ν / (C_crit * B) )
end