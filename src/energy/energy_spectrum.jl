"""
    energy_spectrum(f_p::Vector{<:Real},
                    q::Vector{<:Real},
                    cut::Real,
                    rho::Real,
                    bounds::Vector{<:Real})

Solves the energy integral for each bin. Returns energy per unit mass for each bin.
"""
function energy_spectrum(f_p::Vector{<:Real},
                         q::Vector{<:Real},
                         cut::Real,
                         rho::Real,
                         bounds::Vector{<:Real})

    # we don't need to solve any integrals if there are no CRs
    if iszero(sum(f_p))
        return 0.0
    end

    Nbins = length(f_p)
    CR_E = Vector{Float64}(undef, Nbins)

    # loop over all bins
    for bin = 1:Nbins-1

        # get upper integration boundary
        hbound_proper = bounds[bin+1] > cut ? cut : bounds[bin+1]

        # solve energy integral
        CR_E[bin] = density_integral(bounds[bin], hbound_proper, f_p[bin], q[bin], rho)
    end

    # solve last integral seperately, in case cut > bounds[end]
    hbound_proper = bounds[end] < cut ? cut : bounds[end]

    # solve energy integral
    CR_E[end] = density_integral(bounds[end], hbound_proper, f_p[end], q[end], rho)

    return CR_E
end


"""
    energy_spectrum( f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        rho::Real,
                        par::CRMomentumDistributionConfig)

Solves the energy integral for each bin. Returns energy per unit mass for each bin.

"""
function energy_spectrum(f_p::Vector{<:Real},
                        q::Vector{<:Real},
                        cut::Real,
                        rho::Real,
                        par::CRMomentumDistributionConfig)

    # we don't need to solve any integrals if there are no CRs
    if iszero(sum(f_p))
        return 0.0
    end

    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    # use default computation
    energy_spectrum(f_p, q, cut, rho, bounds)
end

"""
    energy_spectrum( CR::CRMomentumDistribution,
                        rho::Real)

Solves the energy integral for each bin. Returns energy per unit mass for each bin.
"""
function energy_spectrum(CR::CRMomentumDistribution,
                        rho::Real)

    # convert back to primitive variables
    f_p, q, cut = convert(CR)

    # construct boundaries
    bounds = CR.bound[1:2:end]

    # compute energy
    energy_spectrum(f_p, q, cut, rho, bounds)
end