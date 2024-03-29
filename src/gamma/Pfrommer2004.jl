

"""
    α_γ(slope::T) where T

Slope of γ-ray spectrum as a function of the slope of a proton spectrum in momentum space. 
"""
α_γ(slope::T) where T = 4/3 * (slope - 5/2)

"""
    δ_γ(α::T) where T

Shape parameter
"""
δ_γ(α::T) where T = 0.14 * α^(-1.6) + 0.44


"""
    σ_pp(α::T) where T

Scattering cross-section of proton-proton scatterin in ``cm^2``
"""
σ_pp(α::T) where T = 32 * ( 0.96 + exp(4.4 - 2.4α) ) * 1.e-27


"""
    E_γ_powδ(Eγ::T, δ::T) where T

Helper function for the sum in the Integrand.
"""
E_γ_powδ(Eγ::T, δ::T) where {T} = (2Eγ / E_π0)^δ


"""
    E_integrant(Eγ::T, δ::T, α_δ::T) where T

Helper function for the Integrand.
"""
E_integrant(Eγ::T, δ::T, α_δ::T) where T = (E_γ_powδ(Eγ, δ) + E_γ_powδ(Eγ, -δ))^α_δ


"""
    gamma_source_PE04(Eγ::T, f_p::T, q_p::T, nH::T) where {T}

Source function of gamma-ray photons at energy `Eγ` in units of `N_photons erg^-1 s^-1 cm^-3` as given in Pfrommer&Enßlin (2004).

"""
function gamma_source_PE04(Eγ::T, f_p::T, q_p::T, nH::T) where {T}
    
    # slope of γ-spectrum
    α = α_γ(q_p)

    # compute expensive stuff here once
    δ = δ_γ(α)

    # α/δ
    α_δ = -α/δ

    # Werhahn+21, B1
    σ_pp(α_γ(q_p)) * cL * nH * 2^(2 - α) * f_p * 4 / 3α * (1 / E_π0)^α * E_integrant(Eγ, δ, α_δ)

end


"""
    integrate_E_trapez(x_in::Real, θ_steps::Integer=100)

Energy integration step.
"""
function integrate_E_trapez(q::Real, E_steps::Integer=50)

    # slope of γ-spectrum in bin
    α = α_γ(q)

    # prefactor terms
    prefac = σ_pp(α) * 2^(4-α)/3α * mπ_c2^(-α)

    # step size of energy integral between 50GeV and 200GeV
    dE = 150 / E_steps
    E  = 50.0

    # compute expensive stuff here once
    δ = δ_γ(α)

    # α/δ
    α_δ = -α/δ

    # first half step
    ℐ = 0.5 * E_integrant(E, δ, α_δ)
    
    # actual integration
    @inbounds for i ∈ 1:E_steps-1
        E += dE
        ℐ += E_integrant(E, δ, α_δ)
    end

    E += dE
    # last half step
    ℐ += 0.5 * E_integrant(E, δ, α_δ)

    # multiply by step length and prefactor
    ℐ *= dE * prefac

    return ℐ
end



"""
    γ_emissivity_per_bin(norm, slope, bound_low, bound_up, rho)

Computes the γ-ray emissivity per bin. 
Returns number of photons per ``cm^3`` per ``s``.
"""
function γ_emissivity_per_bin(norm, slope, bound_low, bound_up, rho)

    # get CR number density in the bin
    N = density_integral(bound_low, bound_up, norm, slope, rho)

    # get energy integral per bin 
    E = integrate_E_trapez(slope)

    return N * E
end


"""
    γ_emission( f_p::Vector{<:Real},
                q::Vector{<:Real},
                cut::Real,
                rho::Real,
                n_e::Real,
                bounds::Vector{<:Real};
                Emin::Real = 0.78,
                Emax::Real = 1.e3)


"""
function γ_emission(f_p::Vector{<:Real},
                    q::Vector{<:Real},
                    cut::Real,
                    rho::Real,
                    n_e::Real,
                    bounds::Vector{<:Real};
                    pmin::Real = 0.78,
                    pmax::Real = 1.e3)


    # # get minimum and maximum dimensionless momentum
    # # corresponding to energy range
    # pmin = get_p_of_energy(Emin, m_p)
    # pmax = get_p_of_energy(Emax, m_p)

    # if the requested range has already cooled off we can 
    # skip the integrals
    if (pmin > cut) || iszero(sum(f_p))
        return 0.0
    end

    # get the first bin boundary larger than pmin
    start_bin, end_bin = get_start_end_bins(pmin, pmax, cut, bounds)

    # store total emissivity 
    j_γ = 0.0

    # solve first integral in case pmin is within the bin
    if start_bin != 1
        # get upper integration boundary
        hbound_proper = bounds[start_bin] > cut ? cut : bounds[start_bin]

        # construct f_p at pmin
        norm = f_p[start_bin-1] * (bounds[start_bin-1] / pmin)^q[start_bin-1]
        
        # compute emissivity
        j_γ += γ_emissivity_per_bin(norm, q[start_bin-1], pmin, hbound_proper, rho)
        
    end

    # loop over all relevant bins
    @inbounds for bin = start_bin:end_bin-1

        # get upper integration boundary
        hbound_proper = bounds[bin+1] > cut ? cut : bounds[bin+1]

        # solve energy integral
        j_γ += γ_emissivity_per_bin(f_p[bin], q[bin], bounds[bin], hbound_proper, rho)
    end

    # solve last integral seperately, in case cut > bounds[end]
    hbound_proper = bounds[end_bin+1] < cut ? cut : bounds[end_bin+1]

    # check if pmax is smaller than hbound_proper
    hbound_proper = hbound_proper < pmax ? hbound_proper : pmax

    # solve energy integral
    j_γ += γ_emissivity_per_bin(f_p[end_bin], q[end_bin], bounds[end_bin], hbound_proper, rho)

    # add remaining factors and return emissivity
    return j_γ * cL * n_e
end
