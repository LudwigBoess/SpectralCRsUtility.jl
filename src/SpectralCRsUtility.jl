module SpectralCRsUtility

    include("io.jl")
    include("cr_datatypes.jl")
    include("analysis_functions.jl")
    include("get_detailled_data.jl")
    include("synchrotron_kernel.jl")
    include("read_cr_test_data.jl")
    

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
           read_adiabatic_test_data,
           construct_spectrum

           
end # module
