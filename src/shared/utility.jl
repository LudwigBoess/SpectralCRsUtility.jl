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
    find_log_mid(x_start::T, x_end::T) where T

Returns the middle value of two points in log-space.
"""
find_log_mid(x_start::T, x_end::T) where {T} = 10.0^(0.5 * (log10(x_start) + log10(x_end)))


"""
    momentum_bin_boundaries(par::CRMomentumDistributionConfig)

Construct momentum bin boundaries based on the `par` struct.
"""
momentum_bin_boundaries(par::CRMomentumDistributionConfig) = [par.pmin * 10.0^((i - 1) * par.bin_width) for i = 1:par.Nbins+1]

"""
    momentum_bin_centers(par::CRMomentumDistributionConfig)

Construct centers of momentum bins. Disregards spectral cutoff!
"""
momentum_bin_centers(par::CRMomentumDistributionConfig) = [par.pmin * 10.0^((i - 0.5) * par.bin_width) for i = 1:par.Nbins]

"""
    momentum_bin_centers(par::CRMomentumDistributionConfig)

Construct centers of momentum bins. Disregards spectral cutoff!
"""
momentum_bin_centers(bounds::Vector{<:Real}) = [find_log_mid(bounds[i], bounds[i+1]) for i = 1:length(bounds)-1]

"""
    function momentum_bin_centers(bounds::Vector{<:Real}, cut::Real)

Construct centers of momentum bins accounting for spectral cutoff.
"""
function momentum_bin_centers(bounds::Vector{<:Real}, cut::Real)

    Nbins = length(bounds)-1
    bin_centers = Vector{Float64}(undef, Nbins)

    @inbounds for bin = 1:Nbins

        # start of bin
        p_start = bounds[bin]

        # if start of naive bin is above cutoff all other bins are empty
        # -> set to cutoff for nicer plotting
        if p_start > cut
            bin_centers[bin:end] .= cut
            break
        end

        p_end = bounds[bin+1]
        # check if bin is only partially filled
        if p_end > cut
            p_end = cut
        end

        # construct log mid of bin
        bin_centers[bin] = find_log_mid(p_start, p_end)
    end

    return bin_centers
end


"""
    interpolate_spectrum(p, f_p_start, p_start, q)

Interpolates the spectrum to an arbitrary momentum within the bin.
"""
interpolate_spectrum(p, f_p_start, p_start, q) = f_p_start * (p_start / p)^q