"""
    convert(CR::CRMomentumDistribution)

Converts a `CRMomentumDistribution` back into  the primitive variables found in the block.
"""
function convert(CR::CRMomentumDistribution)

    # get number of bins
    Nbins = length(CR.norm) >> 1
    
    # convert back to primitive variables
    f_p = CR.norm[1:2:end]
    q   = [-(log10(CR.norm[i]) - log10(CR.norm[i+1])) / (log10(CR.bound[i]) - log10(CR.bound[i+1]) ) for i = 1:2:2Nbins ]
    cut = CR.bound[end]

    return f_p, q, cut
end
