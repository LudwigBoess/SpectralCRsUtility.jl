"""
    stokes_parameters(f_p::Vector{<:Real},
                      q::Vector{<:Real},
                      cut::Real,
                      B::Vector{<:Real},
                      bounds::Vector{<:Real};
                      ν0::Real = 1.4e9,
                      integrate_pitch_angle::Bool = true)


Computes the Stokes parameters `Q` and `U` (in ``[erg/cm^3/Hz/s]``) for a CR distribution function `f(p)`.

``
j_{\\nu,\\mathrm{pol}} = \\frac{\\sqrt{3} e^3}{c} \\: B_\\parallel \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: \\scrG(x)

Q = j_{\\nu,\\mathrm{pol}} \\sin(2χ)
U = j_{\\nu,\\mathrm{pol}} \\cos(2χ)
``

# Arguments
- `f_p::Vector{<:Real}`:    Spectral Norm for momenta `p`.
- `q::Vector{<:Real}`:      Slopes `q` for momenta `p`.
- `cut::Real`:              Spectral cutoff momentum.
- `B::Vector{<:Real}`:      Magnetic field vector in Gauss.
- `bounds::Vector{<:Real}`: Boundaries of spectral bins 

# Keyword Arguments
- `ν0::Real=1.4e9`:                    Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.

"""
function stokes_parameters( f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            B::Vector{<:Real},
                            bounds::Vector{<:Real};
                            ν0::Real = 1.4e9,
                            integrate_pitch_angle::Bool = true)

    # absolute value of Bfield in image plane
    B_cgs = √(B[1]^2 + B[2]^2)

    # if all norms are 0 -> j_nu = 0!
    if iszero(sum(f_p)) || iszero(B_cgs)
        return 0.0
    end

    # store number of bins 
    Nbins = length(f_p)

    # prefactor to Eq. 17
    # include magnetic field into this
    j_ν_prefac = j_ν_prefac_c * B_cgs

    # if run without pitch angle integration
    # sinθ = 1.0 -> integral factor π/2
    if !integrate_pitch_angle
        j_ν_prefac *= 0.5π
    end

    # storage for synchrotron emissivity
    jν = 0.0

    # find minimum momentum that contributes to synchrotron emission
    p_min_synch = smallest_synch_bright_p(ν0, B_cgs)

    @inbounds for i = 1:Nbins

        # bin integration points
        p_start = bounds[i]
        # check if bin is below cutoff
        if p_start > cut
            break
        end

        p_end = bounds[i+1]
        # check if bin is only partially filled
        if p_end > cut
            p_end = cut
        end

        # if the end of the bin does not contribute to the
        # emission we can skip the bin!
        if p_end < p_min_synch
            continue
        end

        # spectrum integration points
        f_p_start = copy(f_p[i])

        if isnan(f_p_start)
            continue
        end

        # calculate the emissivity of the bin
        jν += emissivity_per_bin(f_p_start, p_start,
                                 p_end, q[i],
                                 B_cgs, ν0, 
                                 𝒢,
                                 integrate_pitch_angle)
    end

    # multiply with constants and Bfield
    jν *= j_ν_prefac
    
    # compute 2χ factors
    sin_2χ = 2.0 * B[1] * B[2] / (B[1] * B[1] + B[2] * B[2])
    cos_2χ = (B[1] * B[1] - B[2] * B[2]) / (B[1] * B[1] + B[2] * B[2])

    # stokes parameters 
    Q = jν * cos_2χ
    U = jν * sin_2χ

    return Q, U
end 

"""
    stokes_parameters(f_p::Vector{<:Real},
                      q::Vector{<:Real},
                      cut::Real,
                      B::Vector{<:Real},
                      par::CRMomentumDistributionConfig;
                      ν0::Real = 1.4e9,
                      integrate_pitch_angle::Bool = fals)

Computes the Stokes parameters `Q` and `U` (in ``[erg/cm^3/Hz/s]``) for a CR distribution function `f(p)`.

``
j_{\\nu,\\mathrm{pol}} = \\frac{\\sqrt{3} e^3}{c} \\: B_\\parallel \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: \\scrG(x)

Q = j_{\\nu,\\mathrm{pol}} \\sin(2χ)
U = j_{\\nu,\\mathrm{pol}} \\cos(2χ)
``

# Arguments
- `f_p::Vector{<:Real}`: Spectral Norm for momenta `p`.
- `q::Vector{<:Real}`:   Slopes `q` for momenta `p`.
- `cut::Real`:           Spectral cutoff momentum.
- `B::Vector{<:Real}`:   Magnetic field vector in Gauss.

# Keyword Arguments
- `ν0::Real=1.4e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
"""
function stokes_parameters( f_p::Vector{<:Real},
                            q::Vector{<:Real},
                            cut::Real,
                            B::Vector{<:Real},
                            par::CRMomentumDistributionConfig;
                            ν0::Real = 1.4e9,
                            integrate_pitch_angle::Bool = true)

    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    # use default computation
    stokes_parameters(f_p, q, cut, B, bounds;
                      ν0, integrate_pitch_angle)
end 

"""
    stokes_parameters( CRe::CRMomentumDistribution,
                           B::Vector{<:Real},
                           par::CRMomentumDistributionConfig;
                           ν0::Real = 1.4e9,
                           integrate_pitch_angle::Bool = false,
                           convert_to_mJy::Bool = false,
                           reduce_spectrum::Bool = true)

Computes the Stokes parameters `Q` and `U` (in ``[erg/cm^3/Hz/s]``) for a `CRMomentumDistribution`.

``
j_{\\nu,\\mathrm{pol}} = \\frac{\\sqrt{3} e^3}{c} \\: B_\\parallel \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: \\scrG(x)

Q = j_{\\nu,\\mathrm{pol}} \\sin(2χ)
U = j_{\\nu,\\mathrm{pol}} \\cos(2χ)
``

# Arguments
- `f_p::Vector{<:Real}`: Spectral Norm for momenta `p`.
- `p::Vector{<:Real}`:   Momenta `p` for number densities.
- `B::Vector{<:Real}`:   Magnetic field vector in Gauss.

# Keyword Arguments
- `ν0::Real=1.4e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``
"""
function stokes_parameters( CRe::CRMomentumDistribution,
                            B::Vector{<:Real},
                            par::CRMomentumDistributionConfig;
                            ν0::Real = 1.4e9,
                            integrate_pitch_angle::Bool = false)

    # absolute value of Bfield in image plane
    B_cgs = √(B[1]^2 + B[2]^2)

    # if all norms are 0 -> j_nu = 0!
    if iszero(sum(CRe.norm)) || iszero(B_cgs)
        return 0.0
    end

    # convert back to primite variables
    f_p, q, cut = convert(CRe)

    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    # use default computation
    stokes_parameters(f_p, q, cut, B, bounds;
                      ν0, integrate_pitch_angle)
end


"""
    stokes_parameters(CRe::CRMomentumDistribution,
                      B::Vector{<:Real},
                      bounds::Vector{<:Real};
                      ν0::Real = 1.4e9,
                      integrate_pitch_angle::Bool = false,
                      convert_to_mJy::Bool = false,
                      reduce_spectrum::Bool = true)

Computes the Stokes parameters `Q` and `U` (in ``[erg/cm^3/Hz/s]``) for a `CRMomentumDistribution`.

``
j_{\\nu,\\mathrm{pol}} = \\frac{\\sqrt{3} e^3}{c} \\: B_\\parallel \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: \\scrG(x)

Q = j_{\\nu,\\mathrm{pol}} \\sin(2χ)
U = j_{\\nu,\\mathrm{pol}} \\cos(2χ)
``

# Arguments
- `f_p::Vector{<:Real}`: Spectral Norm for momenta `p`.
- `p::Vector{<:Real}`:   Momenta `p` for number densities.
- `B::Vector{<:Real}`:   Magnetic field vector in Gauss.

# Keyword Arguments
- `ν0::Real=1.4e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
"""
function stokes_parameters( CRe::CRMomentumDistribution,
                            B::Vector{<:Real},
                            bounds::Vector{<:Real};
                            ν0::Real = 1.4e9,
                            integrate_pitch_angle::Bool = false,
                            convert_to_mJy::Bool = false,
                            reduce_spectrum::Bool = true)

    # absolute value of Bfield in image plane
    B_cgs = √(B[1]^2 + B[2]^2)

    # if all norms are 0 -> j_nu = 0!
    if iszero(sum(CRe.norm)) || iszero(B_cgs)
        return 0.0
    end

    # convert back to primite variables
    f_p, q, cut = convert(CRe)

    # use default computation
    stokes_parameters(f_p, q, cut, B, bounds;
                         ν0, integrate_pitch_angle,
                         convert_to_mJy, reduce_spectrum)
end