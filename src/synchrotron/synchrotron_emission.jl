using SynchrotronKernel

# !!! WRONG !!!
# ! Use SPHtoGrid.jl

"""
    calculate_synch_intensity(CReNorm::Vector{<:Real}, CReSlope::Vector{<:Real}, CReCut::Real,
                                   bounds::Vector{<:Real}, bin_width::Real, 
                                   B::Real, density::Real; 
                                   ν0::Real=1.4e9,
                                   convert_to_mJy::Bool=false,
                                   sum_intensity::Bool=false)

Calculate synchrotron intensity ``I_ν`` in units ``[erg/cm^3/Hz]``.
"""
function calculate_synch_intensity(CReNorm::Vector{<:Real}, CReSlope::Vector{<:Real}, CReCut::Real,
                                   bounds::Vector{<:Real}, bin_width::Real, 
                                   B::Real, density::Real; 
                                   ν0::Real=1.4e9,
                                   energy_density_conversion::Real=0.0, # To do!
                                   convert_to_mJy::Bool=false,
                                   sum_intensity::Bool=false)


    Nbins = size(CReNorm,1)

    F     = Vector{Float64}( undef, Nbins )
    F_mid = Vector{Float64}( undef, Nbins )
    E     = Vector{Float64}( undef, Nbins )
    dE    = Vector{Float64}( undef, Nbins )
    J     = Vector{Float64}( undef, Nbins )
    N     = Vector{Float64}( undef, Nbins )
    N_mid = Vector{Float64}( undef, Nbins )
    x     = Vector{Float64}( undef, Nbins )

    E_min = bounds[1] * m_e * c_l^2
    E_max = bounds[end] * m_e * c_l^2
    di = log(E_max/E_min)/LMB_SPECTRAL_CRs


    @inbounds for i = 1:Nbins

        # calculate energy at current boundary
        E[i]     = bounds[i] * m_e * c_l^2
        x[i]     = ν0 / (nu_c_prefac * E[i]^2 * B)

		# integrate over half a bin
		bound_up = bounds[i] * 10^(bin_width*0.5)
        
        # if the
        if bounds[i] > CReCut
            N[i]     = 0.0
            N_mid[i] = 0.0
        else

            if bound_up > CReCut
                bound_up = CReCut
            end

            N[i]      = energy_integral(bounds[i], bound_up, CReNorm[i],
                                        CReSlope[i], 1.0 ) # set density = 1.0 to get energy density

            N[i]     *= energy_density_conversion

            Norm_mid  = CReNorm[i] * (bound_up/bounds[i])^(-CReSlope[i])
            
            N_mid[i]  = energy_integral(bound_up, bounds[i+1], Norm_mid,
                                        CReSlope[i], 1.0)

            N_mid[i] *= energy_density_conversion
        end
    end

    dE[1] = E[1] - E_min * exp(-1*di)
    for i = 2:Nbins
        dE[i] = E[i] - E[i-1]
    end

    @inbounds for i = 2:Nbins

        K = synchrotron_kernel(x[i])
        F[i] = N[i] * K

        x_mid = 0.5 * (x[i-1] + x[i])
        K_mid = synchrotron_kernel(x_mid)

        F_mid[i] = N_mid[i] * K_mid

        J[i] = dE[i] / 6.0 * (F[i] + F[i-1] + 4*F_mid[i])
    end

    if convert_to_mJy
        J .*= (ν0 * 1.e26)
    end

    if sum_intensity
        return sum(J)
    else
        return J
    end

end