using DelimitedFiles

function read_adiabatic_test_data(fi::String, have_low_cut::Bool=false)

    d = readdlm(fi)

    if have_low_cut
        Nbins = length(d[1:end-3,1])
        pmin = d[1,1]
        pmax = d[1+Nbins, 1]

        low_cut  = d[end-1,1]
        CRpC     = d[end,1]

        CRpS = d[1:Nbins,2]
        CRpN = d[1:Nbins,3]
    else
        Nbins = length(d[1:end-2,1])
        pmin = d[1,1]
        pmax = d[1+Nbins, 1]

        low_cut  = pmin
        CRpC     = d[end,1]

        CRpS = d[1:Nbins,2]
        CRpN = d[1:Nbins,3]
    end

    par = CRMomentumDistributionConfig(pmin, pmax, Nbins, true)

    cr = CRMomentumDistribution(par.Nbins)

    # get zeroth bin
    cr.CRp_bound[1] = low_cut
    cr.CRp_dis[1] = CRpN[1]

    # all other bins
    j = 2
    for i = 1:Nbins-1

        # upper boundary of bin
        cr.CRp_bound[j] = pmin * 10.0^(par.bin_width*i)
        cr.CRp_dis[j] = CRpN[i] * ( cr.CRp_bound[j]/cr.CRp_bound[j-1])^(-CRpS[i])
        if cr.CRp_bound[j] > CRpC/par.mc_p
            cr.CRp_bound[j] = CRpC/par.mc_p
        end

        # lower bound of next bin
        cr.CRp_bound[j+1] = cr.CRp_bound[j]
        cr.CRp_dis[j+1] = CRpN[i+1]
        if cr.CRp_bound[j] == CRpC/par.mc_p
            cr.CRp_bound[j+1] = CRpC/par.mc_p
        end

        j += 2

    end

    # last boundary
    cr.CRp_bound[j] = pmax
    if cr.CRp_bound[j-1] < CRpC/par.mc_p
        cr.CRp_bound[j] = CRpC/par.mc_p
    end
    cr.CRp_dis[j] = CRpN[Nbins] * ( cr.CRp_bound[j]/cr.CRp_bound[j-1])^(-CRpS[Nbins])
    cr.CRp_bound[j+1] = cr.CRp_bound[j]

    return cr
end