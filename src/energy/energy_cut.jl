
"""
    get_p_of_energy(Emin, mass)

Computes the dimensionless momentum corresponding to an energy given in `GeV`for a particle of a given `mass`.
"""
function get_p_of_energy(Emin::Real, mass::Real)
    Emin / ( mass*cL^2 * erg2eV * 1.e-9)
end

"""
    get_p_of_kinetic_energy(Ekin, mass)

Computes the dimensionless momentum of a particle of a given `mass` with kinetic energy `Ekin` in `GeV`,
from `Ekin = (√(p² + 1) - 1) mc²`.
"""
function get_p_of_kinetic_energy(Ekin::Real, mass::Real)
    t = Ekin / ( mass*cL^2 * erg2eV * 1.e-9)
    √(t * (t + 2))
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
                        CR_type::String="e",
                        ultra_relativistic::Bool=false)

Computes the total energy contained in the particle distribution in the energy range between `Emin` and `Emax`.
Energies are given in `GeV`!
Returns energy per unit mass.

The kinetic energy is `T(p) = (√(p² + 1) - 1) mc²`, integrated per bin, and `Emin` and `Emax` are kinetic energies. `ultra_relativistic = true` uses
`T(p) = p mc²` instead.
The spectrum is zero above `cut`. If `cut` lies beyond `bounds[end]` the last bin extends up to `cut`.
"""
function cr_energy_in_range(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            rho::Real,
                            bounds::Vector{<:Real};
                            Emin::Real = 1.0,
                            Emax::Real = 1.e9,
                            CR_type::String="e",
                            cgs=false,
                            ultra_relativistic::Bool=false)

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
    if ultra_relativistic
        pmin = get_p_of_energy(Emin, particle_mass)
        pmax = get_p_of_energy(Emax, particle_mass)
    else
        pmin = get_p_of_kinetic_energy(Emin, particle_mass)
        pmax = get_p_of_kinetic_energy(Emax, particle_mass)
    end

    if cgs
        # the exact kinetic energy needs the momentum in units of mc
        ultra_relativistic || error("cgs = true requires ultra_relativistic = true!")
        bounds = bounds .* (particle_mass * cL)
        cut *= particle_mass * cL
        pmin *= particle_mass * cL
        pmax *= particle_mass * cL
    end

    Nbins = length(f_p)

    # store total energy
    CR_E = 0.0

    for bin = 1:Nbins

        # part of the bin inside the energy range and below the cutoff
        p_low = max(bounds[bin], pmin)
        p_up = bin == Nbins ? cut : min(bounds[bin+1], cut)
        p_up = min(p_up, pmax)

        if p_up <= p_low || iszero(f_p[bin])
            continue
        end

        # construct f_p at the lower integration boundary
        norm = f_p[bin] * (bounds[bin] / p_low)^q[bin]

        # solve energy integral
        if ultra_relativistic
            CR_E += energy_integral(p_low, p_up, norm, q[bin], rho)
        else
            CR_E += energy_integral_dk(bounds[bin], bounds[bin+1], p_low, p_up, norm, q[bin], rho)
        end
    end

    return CR_E

end


"""
    cr_energy_in_range( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        rho::Real,
                        par::CRMomentumDistributionConfig;
                        Emin::Real = 1.0,
                        Emax::Real = 1.e9,
                        CR_type::String="e",
                        ultra_relativistic::Bool=false)

Computes the total energy contained in the particle distribution in the energy range between `Emin` and `Emax`.
Energies are given in `GeV`!
Returns energy per unit mass.
"""
function cr_energy_in_range(f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            rho::Real,
                            par::CRMomentumDistributionConfig;
                            Emin::Real = 1.0,
                            Emax::Real = 1.e9,
                            CR_type::String="e",
                            ultra_relativistic::Bool=false)

    # we don't need to solve any integrals if there are no CRs
    if iszero(sum(f_p))
        return 0.0
    end

    # construct boundaries
    bounds = momentum_bin_boundaries(par)

    # use default computation
    cr_energy_in_range( f_p, q, cut, rho, bounds;
                        Emin, Emax, CR_type, ultra_relativistic)
end

"""
    cr_energy_in_range( CR::CRMomentumDistribution,
                        rho::Real;
                        Emin::Real = 1.0,
                        Emax::Real = 1.e9,
                        CR_type::String="e",
                        ultra_relativistic::Bool=false)

Computes the total energy contained in the particle distribution in the energy range between `Emin` and `Emax`.
Energies are given in `GeV`!
Returns energy per unit mass.
"""
function cr_energy_in_range(CR::CRMomentumDistribution,
                            rho::Real;
                            Emin::Real = 1.0,
                            Emax::Real = 1.e9,
                            CR_type::String="e",
                            ultra_relativistic::Bool=false)

    # convert back to primitive variables
    f_p, q, cut = convert(CR)

    # construct boundaries
    bounds = CR.bound[1:2:end]

    # compute energy
    cr_energy_in_range(f_p, q, cut, rho, bounds; Emin, Emax, CR_type, ultra_relativistic)
end