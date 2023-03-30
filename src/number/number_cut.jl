"""
    cr_number_in_range( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        rho::Real,
                        bounds::Vector{<:Real};
                        Emin::Real = 1.0,
                        Emax::Real = 1.e9,
                        CR_type::String="e")

Computes the total number of CRs contained in the particle distribution in the energy range between `Emin` and `Emax`.
Energies are given in `GeV`!
Returns energy per unit mass.
"""
function cr_number_in_range(f_p::Vector{<:Real},
    q::Vector{<:Real},
    cut::Real,
    rho::Real,
    bounds::Vector{<:Real};
    Emin::Real=1.0,
    Emax::Real=1.e9,
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


    # store total energy 
    CR_N = 0.0

    # solve first integral
    if start_bin != 1
        # get upper integration boundary
        hbound_proper = bounds[start_bin] > cut ? cut : bounds[start_bin]

        # construct f_p at pmin
        norm = f_p[start_bin-1] * (bounds[start_bin-1] / pmin)^q[start_bin-1]

        # solve energy integral
        CR_N = density_integral(pmin, hbound_proper, norm, q[start_bin-1], rho)
    end

    # loop over all relevant bins
    @inbounds for bin = start_bin:end_bin-1

        # get upper integration boundary
        hbound_proper = bounds[bin+1] > cut ? cut : bounds[bin+1]

        # solve energy integral
        CR_N += density_integral(bounds[bin], hbound_proper, f_p[bin], q[bin], rho)
    end

    # solve last integral seperately, in case cut > bounds[end]
    hbound_proper = bounds[end_bin+1] < cut ? cut : bounds[end_bin+1]

    # check if pmax is smaller than hbound_proper
    hbound_proper = hbound_proper < pmax ? hbound_proper : pmax

    # solve energy integral
    CR_N += density_integral(bounds[end_bin], hbound_proper, f_p[end_bin], q[end_bin], rho)

    return CR_N

end


"""
    cr_number_in_range( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        rho::Real,
                        par::CRMomentumDistributionConfig;
                        Emin::Real = 1.0,
                        Emax::Real = 1.e9,
                        CR_type::String="e")

Computes the total number of CRs contained in the particle distribution in the energy range between `Emin` and `Emax`.
Energies are given in `GeV`!
Returns energy per unit mass.
"""
function cr_number_in_range(f_p::Vector{<:Real},
    q::Vector{<:Real},
    cut::Real,
    rho::Real,
    par::CRMomentumDistributionConfig;
    Emin::Real=1.0,
    Emax::Real=1.e9,
    CR_type::String="e")

    # we don't need to solve any integrals if there are no CRs
    if iszero(sum(f_p))
        return 0.0
    end

    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    # use default computation
    cr_number_in_range(f_p, q, cut, rho, bounds;
        Emin, Emax, CR_type)
end

"""
    cr_number_in_range( CR::CRMomentumDistribution,
                        rho::Real;
                        Emin::Real = 1.0,
                        Emax::Real = 1.e9,
                        CR_type::String="e")

Computes the total number of CRs contained in the particle distribution in the energy range between `Emin` and `Emax`.
Energies are given in `GeV`!
Returns energy per unit mass.
"""
function cr_number_in_range(CR::CRMomentumDistribution,
    rho::Real;
    Emin::Real=1.0,
    Emax::Real=1.e9,
    CR_type::String="e")

    # convert back to primitive variables
    f_p, q, cut = convert(CR)

    # construct boundaries
    bounds = CR.bound[1:2:end]

    # compute energy
    cr_number_in_range(f_p, q, cut, rho, bounds; Emin, Emax, CR_type)
end