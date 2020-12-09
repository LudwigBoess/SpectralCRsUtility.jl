using SpectralCRsUtility, Test

@testset "Synchrotron Kernel" begin

    LMB_SPECTRAL_CRs = 24
    CR_pmin     = 1.0
    CR_pmax     = 1.e6
    CR_DSlope   = 1.e-6
    init_slope  = 4.3
    rho0 = 1.0
    bin_width = log10(CR_pmax/CR_pmin) / LMB_SPECTRAL_CRs
    bounds = zeros(LMB_SPECTRAL_CRs+1)
    for i = 1:LMB_SPECTRAL_CRs+1
        bounds[i] = CR_pmin * 10^(bin_width*(i-1))
    end

    CReNorm  = zeros(LMB_SPECTRAL_CRs)
    CReSlope = init_slope .* ones(LMB_SPECTRAL_CRs)

    CReNorm[1] = 1.0
    for i = 2:LMB_SPECTRAL_CRs
        CReNorm[i]  = CReNorm[i-1] * (bounds[i]/bounds[i-1])^(-init_slope)
    end

    @test_nowarn calculate_synch_intensity(CReNorm, CReSlope, bounds, bin_width, 3.0e-6, 1.0)

end