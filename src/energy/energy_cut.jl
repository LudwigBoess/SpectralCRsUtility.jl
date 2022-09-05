
"""
    get_p_of_energy(Emin, mass)

Computes the dimensionless momentum corresponding to an energy given in `GeV`for a particle of a given `mass`.
"""
function get_p_of_energy(Emin, mass)
    Emin / ( mass*c_light^2 * erg2eV * 1.e-9)
end

"""
    get_start_end_bins(pmin, pmax, bounds)
"""
function get_start_end_bins(pmin, pmax, cut, bounds)

    # get total number of bins
    Nbins = length(bounds)-1

    # get the first bin boundary larger than pmin
    start_bin = 1
    @inbounds for i = 1:Nbins
        if bounds[i] < pmin 
            start_bin += 1
        else
            break 
        end
    end

    # get the last bin boundary smaller than pmax
    end_bin = 0
    @inbounds for i = 1:Nbins
        if bounds[i] < pmax && bounds[i] < cut
            end_bin += 1
        end
    end

    start_bin, end_bin
end

"""
    cr_energy_in_range( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        rho::Real,
                        bounds::Vector{<:Real};
                        Emin::Real = 1.0,
                        Emax::Real = 1.e9,
                        CR_type::String="e")

Computes the total energy contained in the particle distribution in the energy range between `Emin` and `Emax`.
Energies are given in `GeV`!
Returns energy in `erg` per unit mass.
"""
function cr_energy_in_range(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            rho::Real,
                            bounds::Vector{<:Real};
                            Emin::Real = 1.0,
                            Emax::Real = 1.e9,
                            CR_type::String="e")

    # we don't need to solve any integrals if there are no CRs
    if iszero(sum(f_p))
        return 0.0
    end

    # select protons or electrons
    if CR_type == "p"
        particle_mass = m_p 
    elseif CR_type == "e"
        particle_mass = m_e
    else
        error("CR_type must be either 'p' or 'e'!")
    end

    # get minimum and maximum dimensionless momentum
    # corresponding to energy range
    pmin = get_p_of_energy(Emin, particle_mass)
    pmax = get_p_of_energy(Emax, particle_mass)

    # if the requested range has already cooled off we can 
    # skip the integrals
    if pmin > cut
        return 0.0
    end

    # get the first bin boundary larger than pmin
    start_bin, end_bin = get_start_end_bins(pmin, pmax, cut, bounds)


    # solve first integral
    if start_bin == 1
        CR_E = 0.0
    else
        # get upper integration boundary
        hbound_proper = bounds[start_bin] > cut ? cut : bounds[start_bin]

        # construct f_p at pmin
        norm = f_p[start_bin-1] * (bounds[start_bin-1] / pmin)^q[start_bin-1]
        
        # solve energy integral
        CR_E = energy_integral(pmin, hbound_proper, norm, q[start_bin-1], rho)
    end
    
    # loop over all relevant bins
    @inbounds for bin = start_bin:end_bin-1

        # get upper integration boundary
        hbound_proper = bounds[bin+1] > cut ? cut : bounds[bin+1]

        # solve energy integral
        CR_E += energy_integral(bounds[bin], hbound_proper, f_p[bin], q[bin], rho)
    end

    # solve last integral seperately, in case cut > bounds[end]
    hbound_proper = bounds[end_bin+1] < cut ? cut : bounds[end_bin+1]

    # check if pmax is smaller than hbound_proper
    hbound_proper = hbound_proper < pmax ? hbound_proper : pmax

    # solve energy integral
    CR_E += energy_integral(bounds[end_bin], hbound_proper, f_p[end_bin], q[end_bin], rho)

    return CR_E

end