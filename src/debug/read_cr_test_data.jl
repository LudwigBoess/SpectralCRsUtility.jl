using DelimitedFiles


function read_test_spectrum(fi::String, have_low_cut::Bool=false)

    d = readdlm(fi)

    if have_low_cut
        Nbins = length(d[1:end-4,1])
        pmin = d[1,1]
        pmax = d[1+Nbins, 1]

        low_cut  = d[end-2,1]
        CRpC     = d[end-1,1]
        
        t        = d[end,1]

        CRpS = d[1:Nbins,2] .|> Float64
        CRpN = d[1:Nbins,3] .|> Float64
    else
        Nbins = length(d[1:end-3,1])
        pmin = d[1,1]
        pmax = d[1+Nbins, 1]

        low_cut  = pmin
        CRpC     = d[end-1,1]

        t        = d[end,1]

        CRpS = d[1:Nbins,2] .|> Float64
        CRpN = d[1:Nbins,3] .|> Float64
    end

    return t, CRMomentumDistribution( CRpN, CRpS, CRpC, pmin, pmax, 1.0, 4 )
end