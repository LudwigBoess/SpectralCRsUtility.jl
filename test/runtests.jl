using SpectralCRsUtility, Test


const ref_cr_norm = [1.663030196617704e21, 8.702622295786366e20, 4.554074542792442e20, 2.3831428493629112e20, 1.247096543231146e20, 6.526045169952111e19, 3.415074597209487e19, 1.787105513290987e19, 9.351907329475914e18, 4.893844826098128e18, 2.5609446648857636e18, 1.3401400758830172e18, 7.012941152587739e17, 3.669865896462585e17, 1.9204375746186806e17, 1.0049632825990325e17, 5.258963960714824e16, 2.7520111847839544e16, 1.4401250165910148e16, 7.536161462127508e15, 3.943666621228124e15, 2.0637172514876385e15, 1.0799414106564088e15, 5.651323841043832e14, 2.957332762796661e14, 1.547569616590072e14, 8.098418102695647e13, 4.237895023461221e13, 2.2176867138903695e13, 1.1605134940197436e13, 6.072956839964599e12, 3.177972937852402e12, 1.663030424859911e12, 8.702623490175813e11, 4.554075167815382e11, 2.3831435035107687e11, 1.2470968855461409e11, 6.526046961283449e10, 3.4150745972094868e10, 1.7871055132909866e10, 9.351907329475916e9, 4.893844826098127e9, 2.5609446648857636e9, 1.3401400758830173e9, 7.01294115258774e8, 3.669865896462585e8, 1.9204375746186805e8, 1.0049632825990324e8, 5.2589639607148245e7, 2.7520111847839545e7, 1.4401250165910147e7, 7.536161462127508e6, 3.9436666212281235e6, 2.0637172514876386e6, 1.0799414106564086e6, 565132.3841043832, 295732.6268748678, 154756.62182605112, 80984.00319267143, 42378.8571740883, 22176.818440413917, 11605.109456320848, 6072.943504274542, 3177.965959297602, 1663.026772988362, 870.2604379964329, 455.40651674586314, 238.31382703329348, 124.70941470288133, 65.26032630646517, 34.15067097999869, 17.871015889612327, 9.351886793493053, 4.893834079636477, 2.560939041271896, 1.3401371330508587, 0.7012925752773109, 0.36698578377532226, 0.1920433357503206, 0.10049610757862379, 0.05258952412483008, 0.02752005141604469, 0.0144012185420054, 0.007536144913366152, 0.003943657961275957, 0.0020637127197423366, 0.001079939039198104, 0.0005651311431223519, 0.0002957326268748678, 0.00015475662182605113, 8.098400319267144e-5, 4.23788571740883e-5, 2.2176818440413917e-5, 1.160510945632085e-5, 6.072943504274542e-6, 3.177965959297602e-6]
const ref_cr_slope = 4.5 .* ones(length(ref_cr_norm))
const ref_cr_cut  = 1.e6
const ref_cr_B    = 5.e-6


@testset "SpectralCRsUtility" begin
    

    pmin = 1.0
    pmax = 1.e6
    par = CRMomentumDistributionConfig(pmin, pmax, length(ref_cr_norm))

    @testset "Spectra" begin
        @testset "CRMomentumDistribution" begin
            CR = CRMomentumDistribution(ref_cr_norm, ref_cr_slope, ref_cr_cut, pmin, pmax, 1.0, 4)
            @test CR.norm[1:2:end] ≈ ref_cr_norm 
        end
    end

    @testset "Synchrotron" begin
        @testset "f_p" begin
            j_ν = synchrotron_emission( ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par, 
                                        ν0 = 1.4e9, integrate_pitch_angle = false, reduce_spectrum = true)

            @test j_ν ≈ 1.2440932091907414e-27

            j_ν = synchrotron_emission( ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par, 
                                        ν0 = 1.4e9, integrate_pitch_angle = true, reduce_spectrum = true)

            @test j_ν ≈ 5.4753081210232675e-28
        end
        @testset "CRMomentumDistribution" begin
            
            norm_spectrum = CRMomentumDistribution(ref_cr_norm, ref_cr_slope, ref_cr_cut, pmin, pmax, 1.0, 4)

            j_ν = synchrotron_emission( norm_spectrum, ref_cr_B, par, 
                                        ν0 = 1.4e9, integrate_pitch_angle = false, reduce_spectrum = true)

            @test j_ν ≈ 1.2440932091907414e-27

            j_ν = synchrotron_emission( norm_spectrum, ref_cr_B, par, 
                                        ν0 = 1.4e9, integrate_pitch_angle = true, reduce_spectrum = true)

            @test j_ν ≈ 5.4753081210232675e-28
        end
    end

    @testset "γ-ray" begin

        @testset "Pfrommer&Enßlin (2004)" begin
            @testset "helper functions" begin
                q = 4.2
                α = SpectralCRsUtility.α_γ(q)
                δ = δ_γ(α)

                @test α ≈ 2.2666666666666666

                @test δ ≈ 0.4778013911548612

                @test SpectralCRsUtility.σ_pp(α) ≈ 4.139054982268097e-23

                @test SpectralCRsUtility.E_γ_powδ(50.0, δ) ≈ 23.505312382037395

                @test SpectralCRsUtility.E_integrant(50.0, δ, -α/δ) ≈ 3.1011953335993916e-7
            end

            @testset "emissivity" begin 

                pmin = 1.0
                pmax = 1.e6
                Nbins = length(ref_cr_norm)
                bin_width = log10(pmax/pmin)/Nbins
                bounds = [pmin * 10.0^((i - 1) * bin_width) for i = 1:Nbins+1]

                @test γ_emission(ref_cr_norm, ref_cr_slope, ref_cr_cut, 
                                1.0, 1.0, bounds) ≈ 564138.2640155746
            end
        end

        @testset "Werhahn et. al. (2021)" begin
            
            @testset "Kafexhio et. al. (2014)" begin

                @testset "cross-sections" begin
                    @testset "σ1π" begin
                        ref_values = [4.939730593781892e-44, 3.349616183683675e-28, 2.9277443259038443e-27, 4.205182345575466e-27, 4.2098467404072576e-27, 4.102902264540402e-27, 4.04135377875269e-27, 4.01908271987257e-27, 4.019322212939915e-27, 4.0305848981430185e-27]
                        Tp = LinRange(Tp_th, 2.0, 10)
                        for i = 1:length(Tp)
                            @test σ1π(Tp[i]) ≈ ref_values[i]
                        end
                    end
                    @testset "σ2π" begin
                        ref_values = Real[0, 0, 5.950065907462581e-30, 3.502162476195644e-29, 2.0111756307349467e-28, 1.0139566623278537e-27, 3.200131373967198e-27, 5.035149203274868e-27, 5.5755592120165185e-27, 5.678577195476546e-27]
                        Tp = LinRange(Tp_th, 2.0, 10)
                        for i = 1:length(Tp)
                            @test σ2π(Tp[i]) ≈ ref_values[i]
                        end
                    end
                end

                @testset "𝓃_π0_K14" begin
                    ref_values = [0.1623951416932937, 0.5441722817877187, 1.1233635727421425, 1.9124386045483552, 2.946840094326904, 4.294155003380654, 5.99847652766958, 8.136370093740226, 10.877581752758365, 14.508457777229449]
                    Tp = 10.0 .^ LinRange(0, 5, 10)
                    for i = 1:length(Tp)
                        @test 𝓃_π0_K14(Tp[i]) ≈ ref_values[i]
                    end
                end
            end

            @testset "Yang et. al. (2018)" begin

                @testet "𝒮" begin
                    @test 𝒮(0.5) ≈ 0.6224593312018546
                end

                @testset "x0" begin
                    @test x0_Y18(1.0) == 0.17
                    @test x0_Y18(3.0) == 0.1
                    @test x0_Y18(6.0) == 0.08
                end

                @testset "α" begin
                    @test α_Y18(0.5) ≈ 0.3667105537514401
                    @test α_Y18(1.2) == 0.36
                    @test α_Y18(7.2) ≈ 0.5624265542500306
                end

                @testset "β" begin
                    @test β_Y18(0.5) ≈ 10.561742124018977
                    @test β_Y18(1.2) == 12.3
                    @test β_Y18(2.5) ≈ 18.500703568578157
                    @test β_Y18(7.2) ≈ 28.5
                end

                @testset "γ" begin
                    @test γ_Y18(0.5) ≈ 1.4342110750288033
                    @test γ_Y18(1.2) == 0.68
                    @test γ_Y18(2.5) ≈ 2.0014101708521315
                end

                @testset "𝒩_low" begin
                    @test 𝒩_low(0.5) ≈ 98.53507611819482
                    @test 𝒩_low(1.2) ≈ 133.3177312831666
                    @test 𝒩_low(2.5) ≈ 301.92303725497464
                    @test 𝒩_low(7.2) ≈ 718.2181602894138
                end

                @testset "𝒩_high" begin
                    @test 𝒩_high(0.5) ≈ 5.308245862802157
                    @test 𝒩_high(1.2) ≈ 3.899779786674582
                    @test 𝒩_high(2.5) ≈ 5.063538597169587
                    @test 𝒩_high(7.2) ≈ 8.479678082265504
                end

                @testset "f_low" begin
                    @test f_low(0.5, 0.5) ≈ 0.25067308096839785
                    @test f_low(0.5, 1.2) ≈ 0.14221547465771092
                    @test f_low(0.5, 2.5) ≈ 0.014504057755778775
                    @test f_low(0.5, 7.2) ≈ 0.00023255732289278092
                end

                @testset "f_high" begin
                    @test f_high(0.5, 0.5) ≈ 1.2956450984182628
                    @test f_high(0.5, 1.2) ≈ 1.387873758732234
                    @test f_high(0.5, 2.5) ≈ 0.930729399594792
                    @test f_high(0.5, 7.2) ≈ 0.5288193022452127
                end

                @testset "f_Y18" begin
                    @test f_high(0.05, 0.5) ≈ 1.065464424851917
                    @test f_high(0.05, 1.2) ≈ 1.2973863168276853
                    @test f_high(0.05, 2.5) ≈ 1.9935155666320428
                    @test f_high(0.05, 7.2) ≈ 4.857609138743655

                    @test f_high(0.5, 0.5) ≈ 0.8205183669119626
                    @test f_high(0.5, 1.2) ≈ 0.8882392055886298
                    @test f_high(0.5, 2.5) ≈ 0.6207630981667199
                    @test f_high(0.5, 7.2) ≈ 0.23139728426253223
                end


            end

            @testset "shared" begin
                @testset "σ_pp_inel" begin
                    @test σ_pp_inel(Tp_th) == 0.0
                    @test σ_pp_inel(1.0)   ≈ 2.2519027368294765e-26
                    @test σ_pp_inel(10.0)  ≈ 2.9469986501180466e-26
                    @test σ_pp_inel(100.0) ≈ 3.1276509767139205e-26
                    @test σ_pp_inel(1.e5)  ≈ 4.7856160923282935e-26
                    @test σ_pp_inel(1.e6)  ≈ 5.719963756550566e-26
                end
            end
        end
    end

end