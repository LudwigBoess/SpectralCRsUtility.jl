"""
    ν_over_ν_crit(p::Real, B::Real, ν::Real)

Computes the fraction ``x = \\frac{ν}{ν_c}`` needed for the synchrotron Kernel.
See Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 19 converted to dimensionless momentum.
"""
ν_over_ν_crit(p::Real, B::Real, ν::Real, sinθ::Real=1.0) = ν / (C_crit_p * B * sinθ * p^2)


"""
    integrate_θ(x_in::Real, synch_F::Function, θ_steps::Integer=50)

Pitch angle integration in Donnert+16, Eq. 17.
`synch_F` can be either the first or second synchrotron function.
"""
function integrate_θ(x_in::Real, synch_F::Function, θ_steps::Integer=50)

    dθ = 0.5π / θ_steps
    θ = 0.0

    # first half step: sin(0) = 0
    K = 0.0

    # actual integration
    @inbounds for i ∈ 1:θ_steps-1
        θ += dθ
        sinθ = sin(θ)
        x = x_in / sinθ
        K += sinθ^2 * synch_F(x)
    end

    # last step: sin(0.5π) = 1
    K += 0.5 * synch_F(x_in)

    # multiply by step length
    K *= dθ

    return K
end


"""
    solve_synchrotron_function( x::Real, synch_F::Function, 
                                integrate_pitch_angle::Bool)

Helper function for solving the pitch angle integration.
"""
function solve_synchrotron_function(x::Real, synch_F::Function, 
                                    integrate_pitch_angle::Bool)
    if integrate_pitch_angle
        return integrate_θ(x, synch_F)
    else
        return synch_F(x)
    end
end

"""
    emissivity_per_bin( f_p_start::Real, p_start::Real,
                        p_end::Real, q::Real,
                        B_cgs::Real, ν0::Real,
                        synch_F::Function,
                        integrate_pitch_angle::Bool)

Calculate the emissivity contained in a spectral bin defined by its start, mid and end values.
"""
function emissivity_per_bin(f_p_start::Real, p_start::Real,
                            p_end::Real, q::Real,
                            B_cgs::Real, ν0::Real,
                            synch_F::Function,
                            integrate_pitch_angle::Bool)

    # x from Eq. 19
    x = ν_over_ν_crit(p_start, B_cgs, ν0)

    # pitch-angle integral
    K = solve_synchrotron_function(x, synch_F, integrate_pitch_angle)

    # energy density at momentum p_start * integrated synchrotron kernel
    F_start = p_start^2 * f_p_start * K

    # middle of bin
    p_mid = find_log_mid(p_start, p_end)
    f_p_mid = interpolate_spectrum(p_mid, f_p_start, p_start, q)

    x = ν_over_ν_crit(p_mid, B_cgs, ν0)
    K = solve_synchrotron_function(x, synch_F, integrate_pitch_angle)

    F_mid = p_mid^2 * f_p_mid * K

    # end of bin 
    x = ν_over_ν_crit(p_end, B_cgs, ν0)
    f_p_end = interpolate_spectrum(p_end, f_p_start, p_start, q)

    K = solve_synchrotron_function(x, synch_F, integrate_pitch_angle)

    F_end = p_end^2 * f_p_end * K

    # bin width
    dp = p_end - p_start

    # store total synchrotron emissivity
    # Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
    return dp / 6 * (F_start + F_end + 4F_mid)
end