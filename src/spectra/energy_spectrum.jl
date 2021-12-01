"""
    getEnergySpectrum( CR_N::Vector{T}, CR_S::Vector{T}, 
                            CR_Cut::T, ρ::T; 
                            par::CRMomentumDistributionConfig, 
                            mode::Integer=3) where T
    
Compute the energy spectrum in ´[10^10 erg/g]´.
"""
function getEnergySpectrum(CR_N::Vector{T}, CR_S::Vector{T},
                            CR_Cut::T, ρ::T;
                            par::CRMomentumDistributionConfig,
                            mode::Integer = 3) where {T<:Real}


    # transform the norm dependent on IO mode
    norm = transform_norm(CR_N, mode)


    bounds = 10.0 .^ collect(log10(par.pmin):par.bin_width:log10(par.pmax))
    bound_up = Vector{Float64}(undef, par.Nbins)

    for i = 1:par.Nbins-1
        bound_up[i] = (bounds[i+1] > CR_Cut) ? CR_Cut : bounds[i+1]
    end

    bound_up[end] = (bounds[end] > CR_Cut) ? CR_Cut : bounds[end]

    bin_centers = Vector{Float64}(undef, par.Nbins)
    energy = Vector{Float64}(undef, par.Nbins)

    @inbounds for i = 1:par.Nbins
        bin_centers[i] = bounds[i] # 10.0^(0.5* (log10(bounds[i]) + log10(bound_up[i])))
        energy[i] = energy_integral(bounds[i], bound_up[i], norm[i], CR_S[i], ρ)

    end

    return CR_EnergySpectrum(bin_centers, energy)
end


"""
    getEnergySpectrum( spec::CRMomentumDistribution, ρ::T ) where T
    
Compute the energy spectrum in ´[10^10 erg/g]´.
"""
function getEnergySpectrum( spec::CRMomentumDistribution, ρ::T ) where T

    Nbins = Int64( 0.5* size(spec.norm,1) )
    bin_centers = Vector{Float64}(undef, Nbins)
    energy      = Vector{Float64}(undef, Nbins)

    k = 1
    @inbounds for i = 1:Nbins
        bin_centers[i] = spec.bound[k] #10.0^(0.5* (log10(spec.bound[k]) + log10(spec.bound[k+1])))
        slope          = log10( spec.bound[k] / spec.bound[k+1] ) /
                         log10( spec.norm[k]  / spec.norm[k+1] )
        energy[i]      = energy_integral(spec.bound[k], spec.bound[k+1], spec.norm[k], slope, ρ)
        k += 2
    end

    return CR_EnergySpectrum(bin_centers, energy)
end
