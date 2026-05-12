using SpectralCRsUtility, Test


const ref_cr_norm = [1.663030196617704e21, 8.702622295786366e20, 4.554074542792442e20, 2.3831428493629112e20, 1.247096543231146e20, 6.526045169952111e19, 3.415074597209487e19, 1.787105513290987e19, 9.351907329475914e18, 4.893844826098128e18, 2.5609446648857636e18, 1.3401400758830172e18, 7.012941152587739e17, 3.669865896462585e17, 1.9204375746186806e17, 1.0049632825990325e17, 5.258963960714824e16, 2.7520111847839544e16, 1.4401250165910148e16, 7.536161462127508e15, 3.943666621228124e15, 2.0637172514876385e15, 1.0799414106564088e15, 5.651323841043832e14, 2.957332762796661e14, 1.547569616590072e14, 8.098418102695647e13, 4.237895023461221e13, 2.2176867138903695e13, 1.1605134940197436e13, 6.072956839964599e12, 3.177972937852402e12, 1.663030424859911e12, 8.702623490175813e11, 4.554075167815382e11, 2.3831435035107687e11, 1.2470968855461409e11, 6.526046961283449e10, 3.4150745972094868e10, 1.7871055132909866e10, 9.351907329475916e9, 4.893844826098127e9, 2.5609446648857636e9, 1.3401400758830173e9, 7.01294115258774e8, 3.669865896462585e8, 1.9204375746186805e8, 1.0049632825990324e8, 5.2589639607148245e7, 2.7520111847839545e7, 1.4401250165910147e7, 7.536161462127508e6, 3.9436666212281235e6, 2.0637172514876386e6, 1.0799414106564086e6, 565132.3841043832, 295732.6268748678, 154756.62182605112, 80984.00319267143, 42378.8571740883, 22176.818440413917, 11605.109456320848, 6072.943504274542, 3177.965959297602, 1663.026772988362, 870.2604379964329, 455.40651674586314, 238.31382703329348, 124.70941470288133, 65.26032630646517, 34.15067097999869, 17.871015889612327, 9.351886793493053, 4.893834079636477, 2.560939041271896, 1.3401371330508587, 0.7012925752773109, 0.36698578377532226, 0.1920433357503206, 0.10049610757862379, 0.05258952412483008, 0.02752005141604469, 0.0144012185420054, 0.007536144913366152, 0.003943657961275957, 0.0020637127197423366, 0.001079939039198104, 0.0005651311431223519, 0.0002957326268748678, 0.00015475662182605113, 8.098400319267144e-5, 4.23788571740883e-5, 2.2176818440413917e-5, 1.160510945632085e-5, 6.072943504274542e-6, 3.177965959297602e-6]
const ref_cr_slope = 4.5 .* ones(length(ref_cr_norm))
const ref_cr_cut  = 1.e6
const ref_cr_B    = 5.e-6
const Bcomponent = √(ref_cr_B^2 / 2)
const ref_cr_Bvec = [Bcomponent, Bcomponent, 0.0]


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

            @test j_ν ≈ 1.0619625420355843e6

            j_ν = synchrotron_emission( ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par, 
                                        ν0 = 1.4e9, integrate_pitch_angle = true, reduce_spectrum = true)

            @test j_ν ≈ 734149.8593302944
        end
        @testset "CRMomentumDistribution" begin
            
            norm_spectrum = CRMomentumDistribution(ref_cr_norm, ref_cr_slope, ref_cr_cut, pmin, pmax, 1.0, 4)

            j_ν = synchrotron_emission( norm_spectrum, ref_cr_B, par, 
                                        ν0 = 1.4e9, integrate_pitch_angle = false, reduce_spectrum = true)

            @test j_ν ≈ 1.0619625420355843e6

            j_ν = synchrotron_emission( norm_spectrum, ref_cr_B, par, 
                                        ν0 = 1.4e9, integrate_pitch_angle = true, reduce_spectrum = true)

            @test j_ν ≈ 734149.8593302943
        end
    end

    @testset "Stokes Parameters" begin
        @testset "f_p" begin
            j_ν = synchrotron_emission(ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par,
                ν0=1.4e9)

            q, u = stokes_parameters(ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_Bvec, par,
                ν0=1.4e9)

            @test iszero(q)
            @test u ≈ -531632.5801020532

            Π = √(q^2 + u^2) / j_ν
            @test Π ≈ 0.7241472206245468

            ψ = 0.5atan(u / q) |> rad2deg
            @test ψ == -45.0
        end
        @testset "CRMomentumDistribution" begin

            norm_spectrum = CRMomentumDistribution(ref_cr_norm, ref_cr_slope, ref_cr_cut, pmin, pmax, 1.0, 4)

            j_ν = synchrotron_emission(norm_spectrum, ref_cr_B, par,
                ν0=1.4e9)

            q, u = stokes_parameters(norm_spectrum, ref_cr_Bvec, par,
                ν0=1.4e9)

            @test iszero(q)
            @test u ≈ -531632.5801020535

            Π = √(q^2 + u^2) / j_ν
            @test Π ≈ 0.7241472206245468

            ψ = 0.5atan(u / q) |> rad2deg
            @test ψ == -45.0
        end
    end

    @testset "γ-ray" begin

        @testset "Pfrommer&Enßlin (2004)" begin
            @testset "helper functions" begin
                q = 4.2
                α = SpectralCRsUtility.α_γ(q)
                δ = SpectralCRsUtility.δ_γ(α)

                @test α ≈ 2.2666666666666666

                @test δ ≈ 0.4778013911548612

                @test SpectralCRsUtility.σ_pp(α) ≈ 4.2030549822680975e-26

                @test SpectralCRsUtility.E_γ_powδ(50.0, δ) ≈ 23.505312382037395

                @test SpectralCRsUtility.E_integrant(50.0, δ, -α/δ) ≈ 3.1011953335993916e-7
            end

            @testset "emissivity" begin 
                @test gamma_source_PE04(1.0, 1.0, 4.5, 1.0) ≈ 3.3709004342715927e-17
            end
        end

        @testset "Werhahn et. al. (2021)" begin
            
            @testset "Kafexhiu et. al. (2014)" begin

                @testset "energy" begin
                    # kinetic energy of proton in GeV
                    @test SpectralCRsUtility.T_p(1.0) ≈ 0.938272088
                    # Ratio between proton kinetic energy and rest energy
                    @test SpectralCRsUtility.θ_p(1.0) ≈ 1.0657889249711967
                    # Squared center of mass energy
                    @test SpectralCRsUtility.𝓈(1.0) ≈ 5.397962220479519
                    # Total energy of the pion in the centre-of-mass system
                    @test SpectralCRsUtility.E_π_CM(1.0) ≈ 0.4077650184727047
                    # Minimum allowed energy for the created pion given in the lab frame
                    @test SpectralCRsUtility.E_π_min(1.0) ≈ 1.0045546908834024
                    # Maximum allowed energy for the created pion given in the lab frame
                    @test SpectralCRsUtility.E_π_max_LAB(2.0) ≈ 1.7130401721061672
                    # Velocity in the center of mass system
                    @test SpectralCRsUtility.β_CM(2.0) ≈ 0.7182781065151783
                    # Pion maximum velocity in the lab frame
                    @test SpectralCRsUtility.β_π_LAB(2.0) ≈ 0.996890937462484
                    # Lorentz factor in the center of mass system 
                    @test SpectralCRsUtility.γ_CM(2.0) ≈ 1.437285262211784
                    # Lorentz factor of the pion in the lab system
                    @test SpectralCRsUtility.γ_π_LAB(2.0) ≈ 12.691358092430386
                    # Minimum γ-ray energy from a pion with kinetic energy `Tp`
                    @test SpectralCRsUtility.E_γ_min(2.0) ≈ 0.0026629745121775988
                    # Maximum γ-ray energy from a pion with kinetic energy `Tp`
                    @test SpectralCRsUtility.E_γ_max(2.0) ≈ 1.7103771975939897
                    # cross-referenced with `naima`
                    @test SpectralCRsUtility.E_π_max_LAB(0.9) ≈ 0.6963603064421705
                    @test SpectralCRsUtility.E_π_max_LAB(2.0) ≈ 1.7130401721061672
                    @test SpectralCRsUtility.E_π_max_LAB(10.0) ≈ 9.586742675463034
                    @test SpectralCRsUtility.E_π_max_LAB(50.0) ≈ 49.543400198671876
                    @test SpectralCRsUtility.E_π_max_LAB(200.0) ≈ 199.5340716232349
                end

                @testset "cross-sections" begin
                    @testset "σ1π" begin
                        ref_values = [4.939730593781892e-44, 3.349616183683675e-28, 2.9277443259038443e-27, 4.205182345575466e-27, 4.2098467404072576e-27, 4.102902264540402e-27, 4.04135377875269e-27, 4.01908271987257e-27, 4.019322212939915e-27, 4.0305848981430185e-27]
                        Tp = LinRange(SpectralCRsUtility.Tp_th, 2.0, 10)
                        for i = 1:length(Tp)
                            @test SpectralCRsUtility.σ1π(Tp[i]) ≈ ref_values[i]
                        end
                    end
                    @testset "σ2π" begin
                        ref_values = Real[0, 0, 5.950065907462581e-30, 3.502162476195644e-29, 2.0111756307349467e-28, 1.0139566623278537e-27, 3.200131373967198e-27, 5.035149203274868e-27, 5.5755592120165185e-27, 5.678577195476546e-27]
                        Tp = LinRange(SpectralCRsUtility.Tp_th, 2.0, 10)
                        for i = 1:length(Tp)
                            @test SpectralCRsUtility.σ2π(Tp[i]) ≈ ref_values[i]
                        end
                    end
                end

                @testset "Fit parameters" begin
                    @testset "α" begin
                        @test SpectralCRsUtility.α_K14(10)  == 1
                        @test SpectralCRsUtility.α_K14(100) == 1/2
                    end
                    @testset "β" begin
                        @test SpectralCRsUtility.β_K14(0.9)    ≈ 3.0771079583046586
                        @test SpectralCRsUtility.β_K14(1.5)    ≈ 2.7423614563386574
                        @test SpectralCRsUtility.β_K14(5.5)    ≈ 4.983149979752324
                        @test SpectralCRsUtility.β_K14(30.5)  == 4.2
                        @test SpectralCRsUtility.β_K14(300.5) == 4.9
                    end
                    @testset "γ" begin
                        @test SpectralCRsUtility.γ_K14(0.5)  == 0
                        @test SpectralCRsUtility.γ_K14(3.5)   ≈ 1.6022183463205815
                        @test SpectralCRsUtility.γ_K14(5.5)   ≈ 1.5220999865015488
                        @test SpectralCRsUtility.γ_K14(25.5) == 1
                    end
                    @testset "κ" begin
                        @test SpectralCRsUtility.κ_K14(10.0) ≈ 3.284251914874226
                    end
                    @testset "μ" begin
                        @test SpectralCRsUtility.μ_K14(10.0) ≈ 0.00013093533160607907
                    end
                    @testset "𝒞" begin
                        @test SpectralCRsUtility.𝒞(2.5) ≈ 0.1848674272287909
                    end   
                    @testset "𝒴_γ" begin
                        @test SpectralCRsUtility.𝒴_γ(2.5) ≈ 2.501821876353361
                    end
                    @testset "𝒳_γ" begin
                        @test SpectralCRsUtility.𝒳_γ(10.0, 2.5) ≈ 0.25041299505090747
                    end
                end

                @testset "Amax" begin
                    @test SpectralCRsUtility.A_max(0.2) == 0
                    @test SpectralCRsUtility.A_max(0.9) ≈ 3.6445288415702625e-26
                    @test SpectralCRsUtility.A_max(1.9) ≈ 6.990647621596554e-26
                    @test SpectralCRsUtility.A_max(5.9) ≈ 1.073022582384389e-25
                end

                @testset "F" begin
                    @test SpectralCRsUtility.F_K14(0.9, 1.0) == 0
                    @test SpectralCRsUtility.F_K14(0.9, 2.0) == 0
                    @test SpectralCRsUtility.F_K14(2.9, 1.0) ≈ 0.04135878590357172
                    @test SpectralCRsUtility.F_K14(2.9, 2.0) ≈ 0.0010219312112906762
                    @test SpectralCRsUtility.F_K14(5.9, 1.0) ≈ 0.07149128657886922
                    @test SpectralCRsUtility.F_K14(5.9, 2.0) ≈ 0.008555466110591579
                end

                @testset "Differential Cross-Section" begin
                    @test SpectralCRsUtility.dσγ_dEγ_K14(2.9, 1.0) ≈ 3.387987951017287e-27
                    @test SpectralCRsUtility.dσγ_dEγ_K14(2.9, 2.0) ≈ 8.37135461058665e-29
                end

                @testset "𝓃_π0_K14" begin
                    ref_values = [0.1623951416932937, 0.5441722817877187, 1.1233635727421425, 1.9124386045483552, 2.946840094326904, 4.294155003380654, 5.99847652766958, 8.136370093740226, 10.877581752758365, 14.508457777229449]
                    Tp = 10.0 .^ LinRange(0, 5, 10)
                    for i = 1:length(Tp)
                        @test SpectralCRsUtility.𝓃_π0_K14(Tp[i]) ≈ ref_values[i]
                    end
                end
            end

            # @testset "Yang et. al. (2018)" begin

            #     @testset "𝒮" begin
            #         @test 𝒮(0.5) ≈ 0.6224593312018546
            #     end

            #     @testset "x0" begin
            #         @test x0_Y18(1.0) == 0.17
            #         @test x0_Y18(3.0) == 0.1
            #         @test x0_Y18(6.0) == 0.08
            #     end

            #     @testset "α" begin
            #         @test α_Y18(0.5) ≈ 0.3667105537514401
            #         @test α_Y18(1.2) == 0.36
            #         @test α_Y18(7.2) ≈ 0.5624265542500306
            #     end

            #     @testset "β" begin
            #         @test β_Y18(0.5) ≈ 10.561742124018977
            #         @test β_Y18(1.2) == 12.3
            #         @test β_Y18(2.5) ≈ 18.500703568578157
            #         @test β_Y18(7.2) ≈ 28.5
            #     end

            #     @testset "γ" begin
            #         @test γ_Y18(0.5) ≈ 1.4342110750288033
            #         @test γ_Y18(1.2) == 0.68
            #         @test γ_Y18(2.5) ≈ 2.0014101708521315
            #     end

            #     @testset "𝒩_low" begin
            #         @test 𝒩_low(0.5) ≈ 98.53507611819482
            #         @test 𝒩_low(1.2) ≈ 133.3177312831666
            #         @test 𝒩_low(2.5) ≈ 301.92303725497464
            #         @test 𝒩_low(7.2) ≈ 718.2181602894138
            #     end

            #     @testset "𝒩_high" begin
            #         @test 𝒩_high(0.5) ≈ 5.308245862802157
            #         @test 𝒩_high(1.2) ≈ 3.899779786674582
            #         @test 𝒩_high(2.5) ≈ 5.063538597169587
            #         @test 𝒩_high(7.2) ≈ 8.479678082265504
            #     end

            #     @testset "f_low" begin
            #         @test f_low(0.5, 0.5) ≈ 0.25067308096839785
            #         @test f_low(0.5, 1.2) ≈ 0.14221547465771092
            #         @test f_low(0.5, 2.5) ≈ 0.014504057755778775
            #         @test f_low(0.5, 7.2) ≈ 0.00023255732289278092
            #     end

            #     @testset "f_high" begin
            #         @test f_high(0.5, 0.5) ≈ 1.2956450984182628
            #         @test f_high(0.5, 1.2) ≈ 1.387873758732234
            #         @test f_high(0.5, 2.5) ≈ 0.930729399594792
            #         @test f_high(0.5, 7.2) ≈ 0.5288193022452127
            #     end

            #     @testset "f_Y18" begin
            #         @test f_high(0.05, 0.5) ≈ 1.065464424851917
            #         @test f_high(0.05, 1.2) ≈ 1.2973863168276853
            #         @test f_high(0.05, 2.5) ≈ 1.9935155666320428
            #         @test f_high(0.05, 7.2) ≈ 4.857609138743655

            #         @test f_high(0.5, 0.5) ≈ 0.8205183669119626
            #         @test f_high(0.5, 1.2) ≈ 0.8882392055886298
            #         @test f_high(0.5, 2.5) ≈ 0.6207630981667199
            #         @test f_high(0.5, 7.2) ≈ 0.23139728426253223
            #     end


            # end

            @testset "shared" begin
                @testset "σ_pp_inel" begin
                    @test SpectralCRsUtility.σ_pp_inel(SpectralCRsUtility.Tp_th) == 0.0
                    @test SpectralCRsUtility.σ_pp_inel(1.0)   ≈ 2.2519027368294765e-26
                    @test SpectralCRsUtility.σ_pp_inel(10.0)  ≈ 2.9469986501180466e-26
                    @test SpectralCRsUtility.σ_pp_inel(100.0) ≈ 3.1276509767139205e-26
                    @test SpectralCRsUtility.σ_pp_inel(1.e5)  ≈ 4.7856160923282935e-26
                    @test SpectralCRsUtility.σ_pp_inel(1.e6)  ≈ 5.719963756550566e-26
                end
            end

            # CR proton spectrum
            norm = [3.3150304033868485e6, 205189.2010208896, 209.70705545666286, 0.171493213402529, 0.00016835520082165043, 1.5618300111348838e-7, 1.4777506904172816e-10, 1.413813951130306e-13]
            slope = [1.7010470628738403, 3.941281318664551, 4.142595291137695, 4.0110182762146, 4.051901340484619, 4.037489414215088, 4.035932540893555, 4.048586368560791]
            cut = 336271.59375
            bounds = [0.1, 0.5623413251903492, 3.1622776601683795, 17.78279410038923, 100.0, 562.3413251903492, 3162.2776601683795, 17782.79410038923, 100000.0]
            nH = 1.25218452056217e-5
            V = 7.086811782958359e69
            d = 3.562040251012591e26

            @testset "Source" begin
                @test gamma_source_pions(norm, slope, cut, bounds, nH, 2.0) ≈ 1.3169307459848811
            end

            @testset "Emissivity" begin
                @test gamma_emissivity_pions(norm, slope, cut, bounds, nH, 2.0) ≈ 2.6338614919697623
            end

            @testset "Luminosity" begin
                @test gamma_luminosity_pions(norm, slope, cut, bounds, nH, V, N_integration_steps=100) ≈ 4.269594666538366e68
            end

            @testset "Flux" begin
                @test gamma_flux_pions(norm, slope, cut, bounds, nH, V, d, N_integration_steps=100) ≈ 1.0130775891406854e17
            end
        end
    end

    @testset "Inverse Compton" begin

        # 1 keV photon energy in erg
        E_X = 1.0e3 / erg2eV
        bounds = momentum_bin_boundaries(par)
        V = 7.086811782958359e69
        d = 3.562040251012591e26

        @testset "CMB seed photons" begin
            # Stefan-Boltzmann u = (4σ/c) T^4 with T_CMB0 = 2.7255 K
            @test u_CMB(0.0) ≈ 4.1749e-13 rtol=1e-4
            # ⟨ε⟩ = 2.7 k_B T at z=0
            @test ε_CMB(0.0) ≈ 1.01604e-15 rtol=1e-4
            # T ∝ (1+z) ⇒ u ∝ (1+z)^4, ε ∝ (1+z)
            @test u_CMB(1.0) / u_CMB(0.0) ≈ 16.0
            @test ε_CMB(1.0) / ε_CMB(0.0) ≈ 2.0
        end

        @testset "Kernel" begin
            # BG70 g(q) = 2q ln q + (1+2q)(1-q)
            @test SpectralCRsUtility.ic_kernel(0.5) ≈ 2*0.5*log(0.5) + (1+1)*(0.5)
            # boundary behaviour
            @test SpectralCRsUtility.ic_kernel(1.0) == 0.0
            @test SpectralCRsUtility.ic_kernel(1.0 + eps()) == 0.0
            @test SpectralCRsUtility.ic_kernel(0.0) == 0.0
            @test SpectralCRsUtility.ic_kernel(-0.1) == 0.0
        end

        @testset "Source" begin
            # TODO: pin to first-run value (replace 0.0 with reported result)
            @test ic_source_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, par, E_X) ≈ 0.0
        end

        @testset "Emissivity" begin
            # By construction: emissivity == E_X * source
            s = ic_source_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, par, E_X)
            @test ic_emissivity_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, bounds, E_X) ≈ E_X * s
        end

        @testset "Luminosity" begin
            # TODO: pin to first-run value
            @test ic_luminosity_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, bounds, V) ≈ 0.0
        end

        @testset "Flux" begin
            # TODO: pin to first-run value
            @test ic_flux_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, bounds, V, d) ≈ 0.0
        end

        @testset "Redshift scaling" begin
            # For Thomson IC by power-law electrons N(γ) ∝ γ^(-s), the source at
            # fixed E_X scales as T^((s+5)/2). Here f_p ∝ p^(-4.5) ⇒ s = 2.5,
            # so q_IC(E_X) ∝ T^3.75 ⇒ z=0→z=1 ratio = 2^3.75 ≈ 13.454.
            # Asymptotic only — needs p_min_ic well inside spectrum at both z.
            s0 = ic_source_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, par, E_X; z=0.0)
            s1 = ic_source_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, par, E_X; z=1.0)
            @test s1 / s0 ≈ 2.0^3.75 rtol=0.03
        end

        @testset "Below threshold" begin
            # Need p_min_ic > cut → no electron can upscatter
            # γ²_min = E_X/(4 ε₀) ⇒ pick E_X so γ_min > ref_cr_cut = 1e6
            ε0 = ε_CMB(0.0)
            E_huge = 4 * (2 * ref_cr_cut)^2 * ε0
            @test ic_source_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, par, E_huge) == 0.0
        end

        @testset "CRMomentumDistribution dispatch" begin
            # Wrapped form must match flat-array form
            CRe = CRMomentumDistribution(ref_cr_norm, ref_cr_slope, ref_cr_cut, pmin, pmax, 1.0, 4)
            s_flat = ic_source_cmb(ref_cr_norm, ref_cr_slope, ref_cr_cut, par, E_X)
            s_wrap = ic_source_cmb(CRe, par, E_X)
            @test s_wrap ≈ s_flat
        end
    end

end