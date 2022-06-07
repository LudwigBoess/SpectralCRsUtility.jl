using SynchrotronKernel

"""
    ν_over_ν_crit(p::Real, B::Real, ν::Real)

Computes the fraction ``x = \\frac{ν}{ν_c}`` needed for the synchrotron Kernel.
See Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 19 converted to dimensionless momentum.
"""
ν_over_ν_crit(p::Real, B::Real, ν::Real, sinθ::Real = 1.0) = ν / (C_crit_p * B * sinθ * p^2)


"""
    integrate_θ(x_in::Real, θ_steps::Integer=50)

Pitch angle integration in Donnert+16, Eq. 17 using Simpson's rule.
"""
function integrate_θ_simpson(x_in::Real, θ_steps::Integer = 50)

    dθ = 0.5π / θ_steps

    F = Vector{Float64}(undef, θ_steps)
    F_mid = Vector{Float64}(undef, θ_steps)

    # sin(0.0) = 0.0
    F[1] = 0.0
    F_mid[1] = sin(0.5dθ)^2 * ℱ(x_in / sin(0.5dθ))
    K = 0.0
    @inbounds for i = 2:θ_steps
        # boundary
        sinθ = sin((i - 1) * dθ)
        x = x_in / sinθ
        F[i] = sinθ^2 * ℱ(x)

        # mid point
        sinθ = sin((i - 0.5) * dθ)
        x = x_in / sinθ
        F_mid[i] = sinθ^2 * ℱ(x)

        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        K += dθ / 6.0 * (F[i] + F[i-1] + 4F_mid[i])
    end
    return K
end


"""
    integrate_θ(x_in::Real, θ_steps::Integer=100)

Pitch angle integration in Donnert+16, Eq. 17.
"""
function integrate_θ(x_in::Real, θ_steps::Integer=128)

    dθ = 0.5π / θ_steps
    K = 0.0
    @inbounds for θ ∈ LinRange(0.0, 0.5π, θ_steps)
        sinθ = sin(θ)
        x  = x_in / sinθ
        K += sinθ^2 * ℱ(x) * dθ
    end
    return K
end


"""
    find_log_mid(x_start::T, x_end::T) where T

Returns the middle value of two points in log-space.
"""
find_log_mid(x_start::T, x_end::T) where T = 10.0^( 0.5 * ( log10(x_start) + log10(x_end) ))


"""
    emissivity( f_p_start::T, p_start::T, 
                f_p_mid::T, p_mid::T, 
                f_p_end::T, p_end::T,
                B_cgs::T, ν0::T, integrate_pitch_angle::Bool) where T

Calculate the emissivity contained in a spectral bin defined by its start, mid and end values.
"""
function emissivity_per_bin(f_p_start::Real, p_start::Real, 
                                     f_p_mid::Real, p_mid::Real, 
                                     f_p_end::Real, p_end::Real,
                                     B_cgs::Real, ν0::Real, integrate_pitch_angle::Bool) where T

    # x from Eq. 19
    x = ν_over_ν_crit(p_start, B_cgs, ν0)

    # pitch-angle integral
    if integrate_pitch_angle
        K = integrate_θ(x)
    else
        K = ℱ(x)
    end

    # energy density at momentum p_start * integrated synchrotron kernel
    F_start = 4π * p_start^2 * f_p_start * K

    # middle of bin
    x = ν_over_ν_crit(p_mid, B_cgs, ν0)

    if integrate_pitch_angle
        K = integrate_θ(x)
    else
        K = ℱ(x)
    end

    F_mid = 4π * p_mid^2  * f_p_mid * K

    # end of bin 
    x = ν_over_ν_crit(p_end, B_cgs, ν0)

    if integrate_pitch_angle
        K = integrate_θ(x)
    else
        K = ℱ(x)
    end

    F_end = 4π * p_end^2 * f_p_end * K

    # bin width
    dp = p_end - p_start

    # store total synchrotron emissivity
    # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
    return dp / 6.0 * (F_start + F_end + 4F_mid)
end



"""
    synchrotron_emission( f_p::Vector{<:Real},
                          q::Vector{<:Real},
                          cut::Real,
                          B_cgs::Real,
                          par::CRMomentumDistributionConfig;
                          ν0::Real = 1.4e9,
                          integrate_pitch_angle::Bool = false,
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
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool = true`:      Return a single value of true, or the spectral components if false

"""
function synchrotron_emission(  f_p::Vector{<:Real},
                                q::Vector{<:Real},
                                cut::Real,
                                B_cgs::Real,
                                par::CRMomentumDistributionConfig;
                                ν0::Real = 1.4e9,
                                integrate_pitch_angle::Bool = false,
                                convert_to_mJy::Bool = false,
                                reduce_spectrum::Bool = true)

    # if all norms are 0 -> j_nu = 0!
    if sum(f_p) == 0.0 || B_cgs == 0.0
        return 0.0
    end

    # prefactor to Eq. 17
    j_ν_prefac = √(3) * qe^3 / c_light

    # include magnetic field into this
    j_ν_prefac *= B_cgs

    # if run without pitch angle integration
    # sinθ = 1.0 -> integral factor π/2
    if !integrate_pitch_angle
        j_ν_prefac *= 0.5π
    end

    # storage array for synchrotron emissivity
    jν = Vector{Float64}(undef, par.Nbins)

    if !reduce_spectrum
        bin_centers = Vector{Float64}(undef, par.Nbins)
    end

    # find minimum momentum that contributes to synchrotron emission
    p_min_synch = smallest_synch_bright_p(ν0, B_cgs)

    @inbounds for i = 1:par.Nbins

        # bin integration points
        p_start = par.pmin * 10.0^((i - 1) * par.bin_width)
        # check if bin is below cutoff
        if p_start > cut
            jν[i:end] .= 0.0
            break
        end

        p_end = par.pmin * 10.0^(i * par.bin_width)
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

        p_mid = find_log_mid(p_start, p_end)

        # spectrum integration points
        f_p_start = f_p[i]
        f_p_mid   = f_p[i] * (p_mid / p_start)^(-q[i])
        f_p_end   = f_p[i] * (p_end / p_start)^(-q[i])


        # calculate the emissivity of the bin
        jν[i] = emissivity_per_bin( f_p_start, p_start, 
                            f_p_mid, p_mid,
                            f_p_end, p_end,
                            B_cgs, ν0, 
                            integrate_pitch_angle)

        if !reduce_spectrum
            bin_centers[i] = p_mid
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
    synchrotron_emission( CRe::CRMomentumDistribution,
                           B_cgs::Real,
                           par::CRMomentumDistributionConfig;
                           ν0::Real = 1.4e9,
                           integrate_pitch_angle::Bool = false,
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
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool = true`:      Return a single value of true, or the spectral components if false

"""
function synchrotron_emission(  CRe::CRMomentumDistribution,
                                B_cgs::Real,
                                par::CRMomentumDistributionConfig;
                                ν0::Real = 1.4e9,
                                integrate_pitch_angle::Bool = false,
                                convert_to_mJy::Bool = false,
                                reduce_spectrum::Bool = true)

    # if all norms are 0 -> j_nu = 0!
    if iszero(sum(CRe.norm)) || iszero(B_cgs)
        return 0.0
    end

    # prefactor to Eq. 17
    j_ν_prefac = √(3) * qe^3 / c_light

    # include magnetic field into this
    j_ν_prefac *= B_cgs

    # if run without pitch angle integration
    # sinθ = 1.0 -> integral factor π/2
    if !integrate_pitch_angle
        j_ν_prefac *= 0.5π
    end

    # storage array for synchrotron emissivity
    jν = Vector{Float64}(undef, par.Nbins)

    if !reduce_spectrum
        bin_centers = Vector{Float64}(undef, par.Nbins)
    end

    # find minimum momentum that contributes to synchrotron emission
    p_min_synch = smallest_synch_bright_p(ν0, B_cgs)
    bin = 0

    @inbounds for i = 1:2:2par.Nbins

        bin += 1

         # if the end of the bin does not contribute to the
        # emission we can skip the bin!
        if CRe.bound[i+1] < p_min_synch
            jν[bin] = 0.0
            continue
        end

        # bin integration points
        p_start = CRe.bound[i]
        p_mid   = find_log_mid(CRe.bound[i], CRe.bound[i+1])
        p_end   = CRe.bound[i+1]

        # spectrum integration points
        f_p_start = CRe.norm[i]
        f_p_mid   = find_log_mid(CRe.norm[i], CRe.norm[i+1])
        f_p_end   = CRe.norm[i+1]

        # calculate the emissivity of the bin
        jν[bin] = emissivity_per_bin( f_p_start, p_start, 
                                      f_p_mid, p_mid,
                                      f_p_end, p_end,
                                      B_cgs, ν0, 
                                      integrate_pitch_angle)

        if !reduce_spectrum
            bin_centers[bin] = p_mid
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
    synchrotron_emission( CRe::CRElectrons,
                          B_cgs::Real,
                          par::CRMomentumDistributionConfig;
                          ν0::Real = 1.4e9,
                          integrate_pitch_angle::Bool = false,
                          convert_to_mJy::Bool = false,
                          reduce_spectrum::Bool = true,
                          CR_norm_factor::Real = 4.428270801560534e21)

Computes the synchrotron emission (in ``[erg/cm^3/Hz/s]``) for a CR distribution function `f(p)` as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

``
j_\\nu(t) = \\frac{\\sqrt{3} e^3}{c} \\: B(t) \\: \\sum\\limits_{i=0}^{N_\\mathrm{bins}} \\:\\int\\limits_0^{\\pi/2} d\\theta  \\text{ sin}^2\\theta \\:  \\int\\limits_{\\hat{p}_\\mathrm{i}}^{\\hat{p}_\\mathrm{i+1}} d\\hat{p} \\:\\: 4\\pi \\hat{p}^2 f(\\hat{p}, t) \\: K(x)
``

# Arguments
- `CRe::CRElectrons`: CR electron spectrum.
- `B_cgs::Real`:         Magnetic field strength (absolute value) in Gauss.

# Keyword Arguments
- `ν0::Real=1.4e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool = true`:      Return a single value of true, or the spectral components if false.
- `CR_norm_factor::Real`:          Unit conversion factor for CR norm

"""
function synchrotron_emission( CRe::CRElectrons,
                               B_cgs::Real,
                               par::CRMomentumDistributionConfig;
                               ν0::Real = 1.4e9,
                               integrate_pitch_angle::Bool = false,
                               convert_to_mJy::Bool = false,
                               reduce_spectrum::Bool = true,
                               CR_norm_factor::Real = 4.428270801560534e21)

    synchrotron_emission(CR_norm_factor .* CRe.Norm, CRe.Slope, CRe.Cut, B_cgs, par;
                         ν0, integrate_pitch_angle, 
                         convert_to_mJy, reduce_spectrum)
end 

