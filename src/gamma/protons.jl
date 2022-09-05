
const pmin_π0 = 1.22 / ( m_p*c_light^2 * erg2eV * 1.e-9)


pmin = 1.22 /( c_light)

function γ_emission(f_p::Vector{<:Real},
                    q::Vector{<:Real},
                    cut::Real)

    hbound_proper = bounds[start_bin] > cut ? cut : bounds[start_bin]
    norm = f_p[start_bin-1] * (bounds[start_bin-1] / pmin)^q[start_bin-1]

end