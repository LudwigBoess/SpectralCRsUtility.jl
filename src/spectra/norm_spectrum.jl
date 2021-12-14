
"""
    norm_spectrum( CR_N::Vector{<:Real}, CR_S::Vector{<:Real}, CR_C::Vector{<:Real}, 
                   pmin::Real, pmax::Real, mc::Real, mode::Integer=3)

Computes the distribution function with double boundaries.
"""
function norm_spectrum(CR_N::Vector{<:Real}, CR_S::Vector{<:Real}, CR_C::Real,
    pmin::Real, pmax::Real, mc::Real, mode::Integer = 3)

    Nbins = size(CR_N, 1)

    # transform the norm dependent on IO mode
    norm = transform_norm(CR_N, mode)

    bin_width = log10(pmax / pmin) / Nbins

    CR_dis = Array{Float64,1}(undef, Int(2 * Nbins))
    CR_bound = Array{Float64,1}(undef, Int(2 * Nbins + 1))

    # get zeroth bin
    CR_bound[1] = pmin
    CR_dis[1] = norm[1]

    # all other bins
    j = 2
    for i = 1:Nbins-1

        # upper boundary of bin
        CR_bound[j] = pmin * 10.0^(bin_width * i)
        CR_dis[j] = norm[i] * (CR_bound[j] / CR_bound[j-1])^(-CR_S[i])
        if CR_bound[j] > CR_C / mc
            CR_bound[j] = CR_C / mc
        end

        # lower bound of next bin
        CR_bound[j+1] = CR_bound[j]
        CR_dis[j+1] = norm[i+1]
        if CR_bound[j] == CR_C / mc
            CR_bound[j+1] = CR_C / mc
        end

        j += 2

    end

    # last boundary
    CR_bound[j] = pmax
    if CR_bound[j-1] < CR_C / mc
        CR_bound[j] = CR_C / mc
    end
    CR_dis[j] = norm[Nbins] * (CR_bound[j] / CR_bound[j-1])^(-CR_S[Nbins])
    CR_bound[j+1] = CR_bound[j]

    return CR_bound, CR_dis

end


function norm_spectrum(CR_N::Vector{T}, CR_S::Vector{T},
                        CR_Cut::T, ρ::T;
                        par::CRMomentumDistributionConfig,
                        mode::Integer = 3) where {T}

    # transform the norm dependent on IO mode
    norm = transform_norm(CR_N, mode)

    bounds = 10.0 .^ collect(log10(par.pmin):par.bin_width:log10(par.pmax))
    bound_up = Vector{Float64}(undef, par.Nbins)

    for i = 1:par.Nbins
        bound_up[i] = (bounds[i+1] > CR_Cut) ? CR_Cut : bounds[i+1]
    end

    bound_up[end] = (bounds[end] > CR_Cut) ? CR_Cut : bounds[end]

    bin_centers = Vector{Float64}(undef, par.Nbins)

    @inbounds for i = 1:par.Nbins
        bin_centers[i] = 10.0^(0.5 * (log10(bounds[i]) + log10(bound_up[i])))
    end

    return CR_NormSpectrum(bin_centers, norm)
end