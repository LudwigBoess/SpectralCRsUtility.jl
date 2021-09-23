module SpectralCRsUtility

    include(joinpath("constants", "constants.jl"))
    include(joinpath("integrals", "number_density.jl"))
    include(joinpath("integrals", "energy.jl"))
    include(joinpath("spectra", "norm_spectrum.jl"))
    include(joinpath("datatypes", "spectra.jl"))
    include(joinpath("io", "io.jl"))
    include(joinpath("spectra", "number_density_spectrum.jl"))
    include(joinpath("spectra", "energy_spectrum.jl"))
    include(joinpath("debug", "get_detailled_data.jl"))
    include(joinpath("debug", "read_cr_test_data.jl"))
    include(joinpath("synchrotron", "synchrotron_emission.jl"))
    include(joinpath("synchrotron", "min_momentum.jl"))
    include(joinpath("plotting", "plot_spectra.jl"))
    

    # datatypes and helper functions for LMB_SPECTRAL_CRs
    export CRShockData,          # datatype to analyse single shocked particle
           readSingleCRShockDataFromOutputFile, # as the name says
           CRMomentumDistributionConfig, # config parameters for momentum distribution function
           CRMomentumDistribution,
           getCRMomentumDistributionFromPartID, # function to get distribution function
           calculateCREnergyInCGS,
           calculateCRNumber,
           get_detailled_shock_data,
           get_detailled_Dpp_data,
           get_detailled_radiative_data,
           get_detailled_adiabatic_data,
           calculate_synch_intensity,
           smallest_synch_bright_p,
           read_adiabatic_test_data,
           construct_spectrum,
           plot

           
end # module
