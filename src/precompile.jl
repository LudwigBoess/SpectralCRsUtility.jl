using SnoopPrecompile    # this is a small dependency

@precompile_setup begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    
    # CR proton spectrum
    p_norm = [3.3150304033868485e6, 205189.2010208896, 209.70705545666286, 0.171493213402529, 0.00016835520082165043, 1.5618300111348838e-7, 1.4777506904172816e-10, 1.413813951130306e-13]
    p_slope = [1.7010470628738403, 3.941281318664551, 4.142595291137695, 4.0110182762146, 4.051901340484619, 4.037489414215088, 4.035932540893555, 4.048586368560791]
    p_cut = 336271.59375
    p_bounds = [0.1, 0.5623413251903492, 3.1622776601683795, 17.78279410038923, 100.0, 562.3413251903492, 3162.2776601683795, 17782.79410038923, 100000.0]
    nH = 1.25218452056217e-5
    V = 7.086811782958359e69
    d = 3.562040251012591e26

    # electrons
    ref_cr_norm = [1.663030196617704e21, 8.702622295786366e20, 4.554074542792442e20, 2.3831428493629112e20, 1.247096543231146e20, 6.526045169952111e19, 3.415074597209487e19, 1.787105513290987e19, 9.351907329475914e18, 4.893844826098128e18, 2.5609446648857636e18, 1.3401400758830172e18, 7.012941152587739e17, 3.669865896462585e17, 1.9204375746186806e17, 1.0049632825990325e17, 5.258963960714824e16, 2.7520111847839544e16, 1.4401250165910148e16, 7.536161462127508e15, 3.943666621228124e15, 2.0637172514876385e15, 1.0799414106564088e15, 5.651323841043832e14, 2.957332762796661e14, 1.547569616590072e14, 8.098418102695647e13, 4.237895023461221e13, 2.2176867138903695e13, 1.1605134940197436e13, 6.072956839964599e12, 3.177972937852402e12, 1.663030424859911e12, 8.702623490175813e11, 4.554075167815382e11, 2.3831435035107687e11, 1.2470968855461409e11, 6.526046961283449e10, 3.4150745972094868e10, 1.7871055132909866e10, 9.351907329475916e9, 4.893844826098127e9, 2.5609446648857636e9, 1.3401400758830173e9, 7.01294115258774e8, 3.669865896462585e8, 1.9204375746186805e8, 1.0049632825990324e8, 5.2589639607148245e7, 2.7520111847839545e7, 1.4401250165910147e7, 7.536161462127508e6, 3.9436666212281235e6, 2.0637172514876386e6, 1.0799414106564086e6, 565132.3841043832, 295732.6268748678, 154756.62182605112, 80984.00319267143, 42378.8571740883, 22176.818440413917, 11605.109456320848, 6072.943504274542, 3177.965959297602, 1663.026772988362, 870.2604379964329, 455.40651674586314, 238.31382703329348, 124.70941470288133, 65.26032630646517, 34.15067097999869, 17.871015889612327, 9.351886793493053, 4.893834079636477, 2.560939041271896, 1.3401371330508587, 0.7012925752773109, 0.36698578377532226, 0.1920433357503206, 0.10049610757862379, 0.05258952412483008, 0.02752005141604469, 0.0144012185420054, 0.007536144913366152, 0.003943657961275957, 0.0020637127197423366, 0.001079939039198104, 0.0005651311431223519, 0.0002957326268748678, 0.00015475662182605113, 8.098400319267144e-5, 4.23788571740883e-5, 2.2176818440413917e-5, 1.160510945632085e-5, 6.072943504274542e-6, 3.177965959297602e-6]
    ref_cr_slope = 4.5 .* ones(length(ref_cr_norm))
    ref_cr_cut = 1.e6
    ref_cr_B = 5.e-6
    pmin = 1.0
    pmax = 1.e6
    par = CRMomentumDistributionConfig(pmin, pmax, length(ref_cr_norm))

    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)

        # synchrotron emission 
        synchrotron_emission(ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par,
            ν0=1.4e9, integrate_pitch_angle=false, reduce_spectrum=true)
        synchrotron_emission(ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par,
            ν0=1.4e9, integrate_pitch_angle=false, reduce_spectrum=false)
        synchrotron_emission(ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par,
            ν0=1.4e9, integrate_pitch_angle=true, reduce_spectrum=true)
        synchrotron_emission(ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, par,
            ν0=1.4e9, integrate_pitch_angle=true, reduce_spectrum=false)

        # pions
        gamma_source_pions(p_norm, p_slope, p_cut, p_bounds, nH, 2.0)
        gamma_emissivity_pions(p_norm, p_slope, p_cut, p_bounds, nH, 2.0)
        gamma_luminosity_pions(p_norm, p_slope, p_cut, p_bounds, nH, V, N_integration_steps=100)
        gamma_flux_pions(p_norm, p_slope, p_cut, p_bounds, nH, V, d, N_integration_steps=100)
    end
end