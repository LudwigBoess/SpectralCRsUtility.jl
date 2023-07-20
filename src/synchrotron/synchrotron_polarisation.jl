"""
    polarisation_fraction(synchrotron_intensity, stokes_Q, stokes_U)
"""
function polarisation_fraction(I_nu::Real, Q::Real, U::Real)
    return √(Q^2 + U^2) / I_nu
end

function polarisation_angle(Q::Real, U::Real; deg::Bool=false)
    
    # definition of polarisation angle
    ψ = 0.5atan(U / Q)
    
    if deg 
        ψ = rad2deg(ψ)
    end

    return ψ
end