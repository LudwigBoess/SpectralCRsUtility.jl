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
