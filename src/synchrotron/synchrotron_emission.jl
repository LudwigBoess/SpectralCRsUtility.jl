"""
    synchrotron_emission( f_p::Vector{<:Real},
                          q::Vector{<:Real},
                          cut::Real,
                          B_cgs::Real,
                          bounds::Vector{<:Real};
                          ν0::Real = 1.4e9,
                          integrate_pitch_angle::Bool = true,
                          convert_to_mJy::Bool = false,
                          reduce_spectrum::Bool = true)

Computes the synchrotron emission (in ``[erg/cm^3/Hz/s]``) for a CR distribution function `f(p)` as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

``
j_\\nu(t) = \\frac{\\sqrt{3} e^3}{c} \\: B(t) \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: K(x)
``

# Arguments
- `f_p::Vector{<:Real}`:    Spectral Norm for momenta `p`.
- `q::Vector{<:Real}`:      Slopes `q` for momenta `p`.
- `cut::Real`:              Spectral cutoff momentum.
- `B_cgs::Real`:            Magnetic field strength (absolute value) in Gauss.
- `bounds::Vector{<:Real}`: Boundaries of spectral bins 

# Keyword Arguments
- `ν0::Real=1.4e9`:                    Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=true`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:        Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool = true`:      Return a single value of true, or the spectral components if false

"""
function synchrotron_emission(  f_p::Vector{<:Real},
                                q::Vector{<:Real},
                                cut::Real,
                                B_cgs::Real,
                                bounds::Vector{<:Real};
                                ν0::Real = 1.4e9,
                                integrate_pitch_angle::Bool = true,
                                convert_to_mJy::Bool = false,
                                reduce_spectrum::Bool = true)

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

    # storage array for synchrotron emissivity
    jν = Vector{Float64}(undef, Nbins)

    if !reduce_spectrum
        bin_centers = Vector{Float64}(undef, Nbins)
    end

    # find minimum momentum that contributes to synchrotron emission
    p_min_synch = smallest_synch_bright_p(ν0, B_cgs)

    @inbounds for i = 1:Nbins

        # bin integration points
        p_start = bounds[i]
        # check if bin is below cutoff
        if p_start > cut
            jν[i:end] .= 0.0
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
            jν[i] = 0.0
            continue
        end

        # spectrum integration points
        f_p_start = copy(f_p[i])

        if isnan(f_p_start)
            jν[i] = 0.0
            continue
        end

        # calculate the emissivity of the bin
        jν[i] = emissivity_per_bin( f_p_start, p_start, 
                                    p_end, q[i],
                                    B_cgs, ν0, 
                                    ℱ,
                                    integrate_pitch_angle)

        if !reduce_spectrum
            bin_centers[i] = find_log_mid(p_start, p_end)
        end
    end


    if convert_to_mJy
        j_ν_prefac *= mJy_factor
    end

    if reduce_spectrum
        return j_ν_prefac * sum(jν)
    else
        return bin_centers, j_ν_prefac .* jν
    end
end 

"""
    synchrotron_emission( f_p::Vector{<:Real},
                          q::Vector{<:Real},
                          cut::Real,
                          B_cgs::Real,
                          par::CRMomentumDistributionConfig;
                          ν0::Real = 1.4e9,
                          integrate_pitch_angle::Bool = true,
                          convert_to_mJy::Bool = false,
                          reduce_spectrum::Bool = true)

Computes the synchrotron emission (in ``[erg/cm^3/Hz/s]``) for a CR distribution function `f(p)` as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

``
j_\\nu(t) = \\frac{\\sqrt{3} e^3}{c} \\: B(t) \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: K(x)
``

# Arguments
- `f_p::Vector{<:Real}`: Spectral Norm for momenta `p`.
- `q::Vector{<:Real}`:   Slopes `q` for momenta `p`.
- `cut::Real`:           Spectral cutoff momentum.
- `B_cgs::Real`:         Magnetic field strength (absolute value) in Gauss.

# Keyword Arguments
- `ν0::Real=1.4e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=true`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool = true`:      Return a single value of true, or the spectral components if false

"""
function synchrotron_emission(  f_p::Vector{<:Real},
                                q::Vector{<:Real},
                                cut::Real,
                                B_cgs::Real,
                                par::CRMomentumDistributionConfig;
                                ν0::Real = 1.4e9,
                                integrate_pitch_angle::Bool = true,
                                convert_to_mJy::Bool = false,
                                reduce_spectrum::Bool = true)

    # if all norms are 0 -> j_nu = 0!
    if iszero(sum(f_p)) || iszero(B_cgs)
        return 0.0
    end

    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    # use default computation
    synchrotron_emission(f_p, q, cut, B_cgs, bounds;
                         ν0, integrate_pitch_angle,
                         convert_to_mJy, reduce_spectrum)
end 

"""
    synchrotron_emission( CRe::CRMomentumDistribution,
                           B_cgs::Real,
                           par::CRMomentumDistributionConfig;
                           ν0::Real = 1.4e9,
                           integrate_pitch_angle::Bool = true,
                           convert_to_mJy::Bool = false,
                           reduce_spectrum::Bool = true)

Computes the synchrotron emission (in ``[erg/cm^3/Hz/s]``) for a `CRMomentumDistribution` as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

``
j_\\nu(t) = \\frac{\\sqrt{3} e^3}{c} \\: B(t) \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: K(x)
``

# Arguments
- `f_p::Vector{<:Real}`: Spectral Norm for momenta `p`.
- `p::Vector{<:Real}`:   Momenta `p` for number densities.
- `B_cgs::Real`:         Magnetic field strength (absolute value) in Gauss.

# Keyword Arguments
- `ν0::Real=1.4e9`:                   Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=true`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool = true`:     Return a single value of true, or the spectral components if false

"""
function synchrotron_emission(  CRe::CRMomentumDistribution,
                                B_cgs::Real,
                                par::CRMomentumDistributionConfig;
                                ν0::Real = 1.4e9,
                                integrate_pitch_angle::Bool = true,
                                convert_to_mJy::Bool = false,
                                reduce_spectrum::Bool = true)

    # if all norms are 0 -> j_nu = 0!
    if iszero(sum(CRe.norm)) || iszero(B_cgs)
        return 0.0
    end

    # convert back to primite variables
    f_p, q, cut = convert(CRe)

    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    # use default computation
    synchrotron_emission(f_p, q, cut, B_cgs, bounds;
                         ν0, integrate_pitch_angle,
                         convert_to_mJy, reduce_spectrum)

end


"""
    synchrotron_emission( CRe::CRMomentumDistribution,
                           B_cgs::Real,
                           bounds::Vector{<:Real};
                           ν0::Real = 1.4e9,
                           integrate_pitch_angle::Bool = true,
                           convert_to_mJy::Bool = false,
                           reduce_spectrum::Bool = true)

Computes the synchrotron emission (in ``[erg/cm^3/Hz/s]``) for a `CRMomentumDistribution` as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

``
j_\\nu(t) = \\frac{\\sqrt{3} e^3}{c} \\: B(t) \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: K(x)
``

# Arguments
- `f_p::Vector{<:Real}`: Spectral Norm for momenta `p`.
- `p::Vector{<:Real}`:   Momenta `p` for number densities.
- `B_cgs::Real`:         Magnetic field strength (absolute value) in Gauss.

# Keyword Arguments
- `ν0::Real=1.4e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=true`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool = true`:      Return a single value of true, or the spectral components if false

"""
function synchrotron_emission(  CRe::CRMomentumDistribution,
                                B_cgs::Real,
                                bounds::Vector{<:Real};
                                ν0::Real = 1.4e9,
    integrate_pitch_angle::Bool=true,
                                convert_to_mJy::Bool = false,
                                reduce_spectrum::Bool = true)

    # if all norms are 0 -> j_nu = 0!
    if iszero(sum(CRe.norm)) || iszero(B_cgs)
        return 0.0
    end

    # convert back to primite variables
    f_p, q, cut = convert(CRe)

    # use default computation
    synchrotron_emission(f_p, q, cut, B_cgs, bounds;
                         ν0, integrate_pitch_angle,
                         convert_to_mJy, reduce_spectrum)

end