function get_num_dens_spectrum(CR_N::Vector{T}, CR_S::Vector{T},
    CR_Cut::T, ρ::T;
    par::CRMomentumDistributionConfig,
    mode::Integer = 3) where {T}

    # transform the norm dependent on IO mode
    norm = transform_norm(CR_N, mode)

    bin_width = log10(par.pmax / par.pmin) / par.Nbins

    bounds = 10.0 .^ collect(log10(par.pmin):par.bin_width:log10(par.pmax))
    bound_up = Vector{Float64}(undef, par.Nbins)

    for i = 1:par.Nbins
        bound_up[i] = (bounds[i+1] > CR_Cut) ? CR_Cut : bounds[i+1]
    end

    bound_up[end] = (bounds[end] > CR_Cut) ? CR_Cut : bounds[end]

    bin_centers = Vector{Float64}(undef, par.Nbins)
    density = Vector{Float64}(undef, par.Nbins)

    @inbounds for i = 1:par.Nbins
        bin_centers[i] = 10.0^(0.5 * (log10(bounds[i]) + log10(bound_up[i])))
        density[i] = density_integral(bounds[i], bound_up[i], norm[i], CR_S[i], ρ)
    end

    return CR_NumSpectrum(bin_centers, density)
end