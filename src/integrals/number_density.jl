@inline function density_integral(bound_low::Real, bound_up::Real,
                                  norm::Real, slope::Real, ρ::Real)

    # density integral (eq. 9 M01)
    nb = 4π * norm * bound_low^3 / ρ
    density = nb * ( (bound_up/bound_low)^(3 - slope ) - 1 ) / ( 3 - slope )

    if ( 3 - slope_soft ) < slope < ( 3 + slope_soft )
        slope_var = (slope - 3)/slope_soft
        density2 = nb * log10(bound_up/bound_low)
        if slope_var != 0
            density = density * slope_var + density2 * ( 1 - slope_var )
        else
            density = density2
        end
    end

    return density
end



