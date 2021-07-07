using SynchrotronKernel


function calculate_synch_intensity(CReNorm::Vector{<:Real}, CReSlope::Vector{<:Real}, CReCut::Real,
                                   bounds::Vector{<:Real}, bin_width::Real, 
                                   B::Real, density::Real; 
                                   ν0::Real=1.4e9,
                                   convert_to_mJy::Bool=false,
                                   sum_intensity::Bool=false)

    m_e = 9.10953e-28
    c_l = 2.9979e10
    qe  = 4.8032e-10
    #erg2eV = 6.242e+11

    #E_cntr_prefac = m_e*c_l*c_l * sqrt(2.0/3.0 * (2π * m_e*c_l)/qe /0.29)
    nu_c_prefac = 3.0 * qe / (4π * m_e^3 * c_l^5)
    #j_nu_prefac = qe * qe * qe * sqrt(3.0) / (m_e * c_l * c_l)

    LMB_SPECTRAL_CRs = size(CReNorm,1)

    # F     = Vector{Float64}(undef, LMB_SPECTRAL_CRs)
    # F_mid = Vector{Float64}(undef, LMB_SPECTRAL_CRs)
    # E     = Vector{Float64}(undef, LMB_SPECTRAL_CRs)
    # dE    = Vector{Float64}(undef, LMB_SPECTRAL_CRs)
    # J     = Vector{Float64}(undef, LMB_SPECTRAL_CRs)
    # N     = Vector{Float64}(undef, LMB_SPECTRAL_CRs)
    # N_mid = Vector{Float64}(undef, LMB_SPECTRAL_CRs)
    # x     = Vector{Float64}(undef, LMB_SPECTRAL_CRs)

    F     = zeros(LMB_SPECTRAL_CRs)
    F_mid = zeros(LMB_SPECTRAL_CRs)
    E     = zeros(LMB_SPECTRAL_CRs)
    dE    = zeros(LMB_SPECTRAL_CRs)
    J     = zeros(LMB_SPECTRAL_CRs)
    N     = zeros(LMB_SPECTRAL_CRs)
    N_mid = zeros(LMB_SPECTRAL_CRs)
    x     = zeros(LMB_SPECTRAL_CRs)

    E_min = bounds[1] * m_e * c_l^2
    E_max = bounds[end] * m_e * c_l^2
    di = log(E_max/E_min)/LMB_SPECTRAL_CRs


    for i = 1:LMB_SPECTRAL_CRs
        E[i]     = bounds[i] * m_e * c_l^2
        x[i]     = ν0 / (nu_c_prefac * E[i]^2 * B)

		# old version
        # N[i]     = CReNorm[i]
        # bound_up = bounds[i] * 10^(bin_width*0.5)
        # N_mid[i] = CReNorm[i] * (bound_up/bounds[i])^(-CReSlope[i])

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

            N[i]     = density_integral(bounds[i], bound_up, CReNorm[i],
                                        CReSlope[i], density)

            Norm_mid = CReNorm[i] * (bound_up/bounds[i])^(-CReSlope[i])
            N_mid[i] = density_integral(bound_up, bounds[i+1], Norm_mid,
                                        CReSlope[i], density)
        end
    end

    dE[1] = E[1] - E_min * exp(-1*di)
    for i = 2:LMB_SPECTRAL_CRs
        dE[i] = E[i] - E[i-1]
    end

    for i = 2:LMB_SPECTRAL_CRs

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
