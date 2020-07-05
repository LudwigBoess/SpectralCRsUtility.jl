using SpectralCRsUtility, Test

@testset "Synchrotron Kernel" begin

    @test synchrotron_kernel(1.e-8) ≈ 0.0046310000727714734
    @test synchrotron_kernel(  1.0) ≈ 0.6514228153553652
    @test synchrotron_kernel(  3.5) ≈ 0.08268719536694746
    @test synchrotron_kernel(  4.5) ≈ 0.03357441971502366
    @test synchrotron_kernel( 10.0) ≈ 0.00019223826430086885
    @test synchrotron_kernel( 50.0) ≈ 1.734785203976593e-21
    @test synchrotron_kernel(100.0) ≈ 4.697593665922202e-43
    @test synchrotron_kernel( 1.e6) == 0.0

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