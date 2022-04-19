"""   
    struct CRMomentumDistributionConfig{T1,T2}
        pmin::T1
        pmax::T1
        Nbins::T2
        bin_width::T1
        mp::T1
        me::T1
        c::T1
        mc_e::T1
        mc_p::T1
    end
    
Config object to obtain distribution spectrum in momentum space.
"""
struct CRMomentumDistributionConfig
    pmin::Float64
    pmax::Float64
    Nbins::Int64
    bin_width::Float64
    mp::Float64
    me::Float64
    c::Float64
    mc_e::Float64
    mc_p::Float64

    function CRMomentumDistributionConfig(pmin::Real = 0.0, pmax::Real = 0.0,
        Nbins::Integer = 24, mode::Integer = 3)

        # the original version used cgs units
        if mode == 1
            CNST_ME = 9.1095e-28
            CNST_MP = 1.6726e-24
            CNST_C = 2.9979e10
        else
            CNST_ME = 1.0
            CNST_MP = 1.0
            CNST_C = 1.0
        end

        MCe = CNST_ME * CNST_C
        MCp = CNST_MP * CNST_C
        bin_width = log10(pmax / pmin) / Nbins

        new(pmin, pmax, Nbins, bin_width, CNST_MP, CNST_ME, CNST_C, MCe, MCp)

    end
end


"""
    struct CRMomentumDistribution
        bound::Vector{Float64}
        norm::Vector{Float64}
    end

Distribution function reconstruction.
"""
struct CRMomentumDistribution
    bound::Vector{Float64}
    norm::Vector{Float64}
end

"""
    CRMomentumDistribution(CR_N::Vector{<:Real}, CR_S::Vector{<:Real}, CR_C::Real,
                                    pmin::Real, pmax::Real, mc_p::Real)

Function to construct the `CRMomentumDistribution` via properties from a snapshot file.
"""
function CRMomentumDistribution(CR_N::Vector{<:Real}, CR_S::Vector{<:Real}, CR_C::Real,
    pmin::Real, pmax::Real, mc_p::Real, mode::Integer = 3)

    bound, norm = norm_spectrum(CR_N, CR_S, CR_C, pmin, pmax, mc_p, mode)
    CRMomentumDistribution(bound, norm)
end

# """
#     CRMomentumDistribution( CR::AbstractCRSpectrum,
#                             pmin::Real, pmax::Real )

# Function to construct the `CRMomentumDistribution` via an `AbstractCRSpectrum` from SpectralFkpSolver.jl.
# """
# function CRMomentumDistribution( CR::AbstractCRSpectrum,
#                                  pmin::Real, pmax::Real )

#     bound, norm = norm_spectrum(CR.Norm, CR.Slope, CR.Cut, pmin, pmax, 1.0, 4)
#     CRMomentumDistribution(bound, norm)
# end


# """
#     CRMomentumDistribution( CR::AbstractCRSpectrum,
#                             par::RunParameters )

# Function to construct the `CRMomentumDistribution` via an `AbstractCRSpectrum` from SpectralFkpSolver.jl.
# """
# function CRMomentumDistribution( CR::AbstractCRSpectrum,
#                                  par::RunParameters )

#     bound, norm = norm_spectrum(CR.Norm, CR.Slope, CR.Cut, par.bounds[1], par.bounds[end], 1.0, 4)
#     CRMomentumDistribution(bound, norm)
# end


"""
    struct CR_NormSpectrum
        bin::Vector{Float64}
        energy::Vector{Float64}
    end

Stores the norm spectrum.
"""
struct CR_NormSpectrum
    bin::Vector{Float64}
    norm::Vector{Float64}
end


"""
    struct CR_EnergySpectrum
        bin::Vector{Float64}
        energy::Vector{Float64}
    end

Stores the energy spectrum.
"""
struct CR_EnergySpectrum
    bin::Vector{Float64}
    energy::Vector{Float64}
end


"""
    struct CR_NumSpectrum
        bin::Vector{Float64}
        density::Vector{Float64}
    end

Stores the number density spectrum.
"""
struct CR_NumSpectrum
    bin::Vector{Float64}
    density::Vector{Float64}
end