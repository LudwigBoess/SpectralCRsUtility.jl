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

"""
    get_bounds(par::CRMomentumDistributionConfig)

Construct momentum bin boundaries based on the `par` struct.
"""
function get_bounds(par::CRMomentumDistributionConfig)
    [par.pmin * 10.0^((i - 1) * par.bin_width) for i = 1:par.Nbins+1]
end

"""
    find_log_mid(x_start::T, x_end::T) where T

Returns the middle value of two points in log-space.
"""
find_log_mid(x_start::T, x_end::T) where T = 10.0^( 0.5 * ( log10(x_start) + log10(x_end) ))


"""
    interpolate_spectrum(p, f_p_start, p_start, q)

Interpolates the spectrum to an arbitrary momentum within the bin.
"""
interpolate_spectrum(p, f_p_start, p_start, q) = f_p_start * (p_start / p)^q