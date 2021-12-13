function transform_norm(norm_in, mode)

    Nbins = size(norm_in, 1)
    norm = Vector{Float64}(undef, Nbins)

    @inbounds for i = 1:Nbins
        if mode == 1
            norm[i] = norm_in[i] * 1.e20
        elseif mode == 2
            norm[i] = norm_in[i] * 1.e-20
        elseif mode == 3
            norm[i] = 10.0^norm_in[i]
        elseif mode == 4
            norm[i] = norm_in[i]
        end
    end

    return norm

end

function construct_bin_centers(par, cut)

    bin = Vector{Float64}(undef, par.Nbins)

    @inbounds for i = 1:par.Nbins
    
        # beginning of bin
        p_start = par.pmin * 10.0^((i - 1) * par.bin_width)
    
        p_end = par.pmin * 10.0^(i * par.bin_width)
        if p_end > cut
            p_end = cut
        end

        bin[i] = 10.0^( 0.5 * ( log10(p_start) + log10(p_end) ) )
    end

    return bin
end