using SynchrotronKernel

"""
    Constants for synchrotron
"""

const c_light    = 2.9979e10
const m_e        = 9.10953e-28
const q_e        = 1.602176487e-20 * c_light
const C_crit     = 3q_e / ( 4π * m_e * c_light ) # Donnert+16, MNRAS 462, 2014–2032 (2016), Eg. 20 
                                                        #  -> converted to dimensionless momentum
const γ          = 5.0/3.0
const mJy_factor = 1.e26       # conversion factor from [erg/cm^3/Hz/s] to mJy/cm.


"""
    ν_over_ν_crit(p::Real, B::Real, ν::Real)

Computes the fraction ``x = \frac{ν}{ν_c}`` needed for the synchrotron Kernel.
See Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 19.
"""
ν_over_ν_crit(p::Real, B::Real, ν::Real, sinθ::Real = 1.0) = ν / (C_crit * B * sinθ * p^2)



"""
    integrate_θ(x_in::Real, θ_steps::Integer=100)

Pitch angle integration in Donnert+16, Eq. 17.
"""
function integrate_θ(x_in::Real, θ_steps::Integer = 64)

    dθ = 0.5π / θ_steps
    K = 0.0
    @inbounds for θ ∈ LinRange(0.0, 0.5π, θ_steps)
        sinθ = sin(θ)
        x = x_in / sinθ
        K += sinθ^2 * synchrotron_kernel(x) * dθ
    end
    return K
end

"""
    integrate_θ(x_in::Real, θ_steps::Integer=100)

Pitch angle integration in Donnert+16, Eq. 17 using Simpson's rule.
"""
function integrate_θ_simpson(x_in::Real, θ_steps::Integer = 50)

    dθ = 0.5π / θ_steps
    #θ = LinRange(0.0, 0.5π, θ_steps)

    F = Vector{Float64}(undef, θ_steps)
    F_mid = Vector{Float64}(undef, θ_steps)

    # sin(0.0) = 0.0
    F[1] = 0.0
    F_mid[1] = sin(0.5dθ)^2 * synchrotron_kernel(x_in / sin(0.5dθ))
    K = 0.0
    @inbounds for i = 2:θ_steps
        # boundary
        sinθ = sin((i - 1) * dθ)
        x = x_in / sinθ
        F[i] = sinθ^2 * synchrotron_kernel(x)

        # mid point
        sinθ = sin((i - 0.5) * dθ)
        x = x_in / sinθ
        F_mid[i] = sinθ^2 * synchrotron_kernel(x)

        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        K += dθ / 6.0 * (F[i] + F[i-1] + 4F_mid[i])
    end
    return K
end


"""
    synchrotron_emission( f_p::Vector{<:Real}, 
                                   q::Vector{<:Real},
                                   cut::Real,
                                   B_cgs::Real;
                                   ν0::Real=1.44e9
                                   integrate_pitch_angle::Bool=false,
                                   convert_to_mJy::Bool=false,
                                   reduce_spectrum::Bool = true)

Computes the synchrotron emission (in ``[erg/cm^3/Hz/s]``) for a CR spectrum as described in Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 17.

# Arguments
- `f_p::Vector{<:Real}`: Spectral Norm of bins in cgs units.
- `q::Vector{<:Real}`:   Slopes of bins.
- `cut::Real`:           Spectral cutoff.
- `B_cgs::Real`:         Magnetic field strength (absolute value) in Gauss.
- `par::CRMomentumDistributionConfig`: Parameters of spectrum

# Keyword Arguments
- `ν0::Real=1.44e9`:                  Observation frequency in ``Hz``.
- `integrate_pitch_angle::Bool=false`: Explicitly integrates over the pitch angle. If `false` assumes ``sin(θ) = 1``.
- `convert_to_mJy::Bool=false`:       Convert the result from ``[erg/cm^3/Hz/s]`` to ``mJy/cm``.
- `reduce_spectrum::Bool=true`:       Return result as spectrum if set to false, otherwise calculate total emission.
"""
function synchrotron_emission( f_p::Vector{<:Real},
                               q::Vector{<:Real},
                               cut::Real,
                               B_cgs::Real,
                               par::CRMomentumDistributionConfig;
                               ν0::Real = 1.44e9,
                               integrate_pitch_angle::Bool = false,
                               convert_to_mJy::Bool = false,
                               reduce_spectrum::Bool = true )

    # prefactor to Eq. 17
    j_ν_prefac = q_e * √(3) / (m_e * c_light^2)

    # include magnetic field into this
    j_ν_prefac *= B_cgs

    # if run without pitch angle integration
    # sinθ = 1.0 -> integral factor π/2
    if !integrate_pitch_angle
        j_ν_prefac *= 0.5π
    end

    # storage array for synchrotron emissivity
    jν = Vector{Float64}(undef, par.Nbins)

    # loop over spectrum
    @inbounds for i = 1:par.Nbins
    
        # beginning of bin
        p_start = par.pmin * 10.0^((i - 1) * par.bin_width)
    
        # check if bin is below cutoff
        if p_start > cut
            # set remaining bins to zero
            jν[i:end] = 0
            break
        end
    
        # x from Eq. 19
        x = ν_over_ν_crit(p_start, B_cgs, ν0)
        # pitch-angle integral
        if integrate_pitch_angle
            K = integrate_θ_simpson(x)
        else
            K = synchrotron_kernel(x)
        end
        
        # energy density at momentum p * integrated synchrotron kernel
        F_start = 4π * p_start^2 * f_p[i] * K
    
        # middle of bin
        p_mid = par.pmin * 10.0^((i - 1 / 2) * par.bin_width)
        f_p_mid = f_p[i] * (p_mid / p_start)^(-q[i])
    
        # x from Eq. 19
        x = ν_over_ν_crit(p_mid, B_cgs, ν0)
    
        if integrate_pitch_angle
            K = integrate_θ_simpson(x)
        else
            K = synchrotron_kernel(x)
        end
        
        F_mid = 4π * p_mid^2 * f_p_mid * K
    
        # end of bin 
        p_end = par.pmin * 10.0^(i * par.bin_width)
        if p_end > cut
            p_end = cut
        end
    
        # construct norm at end of spectrum
        f_p_end = f_p[i] * (p_end / p_start)^(-q[i])
    
        # x from Eq. 19
        x = ν_over_ν_crit(p_end, B_cgs, ν0)
    
        if integrate_pitch_angle
            K = integrate_θ_simpson(x)
        else
            K = synchrotron_kernel(x)
        end
    
        F_end = 4π * p_end^2 * f_p_end * K
    
        # bin width
        dp = p_end - p_start
    
        # store total synchrotron emissivity
        # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
        jν[i] = dp / 6.0 * (F_start + F_end + 4F_mid)
    end


    if convert_to_mJy
        j_ν_prefac *= mJy_factor
    end

    if reduce_spectrum
        # return total emissivity
        return j_ν_prefac * sum(jν)
    else
        # return emissivity per bin
        return j_ν_prefac .* jν
    end

end 
