
"""
    norm_spectrum( CR_N::Vector{<:Real}, CR_S::Vector{<:Real}, CR_C::Vector{<:Real}, 
                   pmin::Real, pmax::Real, mc::Real, mode::Integer=3)

Computes the distribution function with double boundaries.
"""
function norm_spectrum( CR_N::Vector{<:Real}, CR_S::Vector{<:Real}, CR_C::Real, 
                        pmin::Real, pmax::Real, mc::Real, mode::Integer=3)
    
    Nbins = size(CR_N,1)

    norm = Vector{Float64}(undef, Nbins)
    
    @inbounds for i = 1:Nbins
        if mode == 1
            norm[i] = CR_N[i] * 1.e20
        elseif mode == 2
            norm[i] = CR_N[i] * 1.e-20
        elseif mode == 3
            norm[i] = 10.0^CR_N[i]
        elseif mode == 4
            norm[i] = CR_N[i]
        end
    end
    
    bin_width = log10(pmax/pmin) / Nbins

    CR_dis   = Array{Float64,1}(undef, Int(2*Nbins))
    CR_bound = Array{Float64,1}(undef, Int(2*Nbins + 1))

    # get zeroth bin
    CR_bound[1] = pmin
    CR_dis[1] = norm[1]

    # all other bins
    j = 2
    for i = 1:Nbins-1

        # upper boundary of bin
        CR_bound[j] = pmin * 10.0^(bin_width*i)
        CR_dis[j] = norm[i] * ( CR_bound[j]/CR_bound[j-1])^(-CR_S[i])
        if CR_bound[j] > CR_C/mc
            CR_bound[j] = CR_C/mc
        end

        # lower bound of next bin
        CR_bound[j+1] = CR_bound[j]
        CR_dis[j+1] = norm[i+1]
        if CR_bound[j] == CR_C/mc
            CR_bound[j+1] = CR_C/mc
        end

        j += 2

    end

    # last boundary
    CR_bound[j] = pmax
    if CR_bound[j-1] < CR_C/mc
        CR_bound[j] = CR_C/mc
    end
    CR_dis[j]     = norm[Nbins] * ( CR_bound[j]/CR_bound[j-1])^(-CR_S[Nbins])
    CR_bound[j+1] = CR_bound[j]

    return CR_bound, CR_dis

end