module SpectralCRsUtility

using SynchrotronKernel

include(joinpath("shared", "constants.jl"))
include(joinpath("integrals", "number_density.jl"))
include(joinpath("integrals", "energy.jl"))
include(joinpath("datatypes", "spectra.jl"))
include(joinpath("spectra", "norm_spectrum.jl"))
include(joinpath("shared", "utility.jl"))
include(joinpath("io", "io.jl"))
include(joinpath("spectra", "utility.jl"))
include(joinpath("spectra", "number_density_spectrum.jl"))
include(joinpath("spectra", "energy_spectrum.jl"))
include(joinpath("debug", "get_detailled_data.jl"))
include(joinpath("debug", "read_cr_test_data.jl"))
include(joinpath("synchrotron", "constants.jl"))
include(joinpath("synchrotron", "min_momentum.jl"))
include(joinpath("synchrotron", "shared_functions.jl"))
include(joinpath("synchrotron", "synchrotron_emission.jl"))
include(joinpath("synchrotron", "synchrotron_polarisation.jl"))
include(joinpath("synchrotron", "stokes_parameters.jl"))
include(joinpath("energy", "energy_cut.jl"))
include(joinpath("energy", "energy_spectrum.jl"))
include(joinpath("number", "number_cut.jl"))
include(joinpath("gamma", "Pfrommer2004.jl"))
include(joinpath("gamma", "energy_computation.jl"))
include(joinpath("gamma", "cross_sections.jl"))
include(joinpath("gamma", "Kafexhiu2014.jl"))
include(joinpath("gamma", "Yang2018.jl"))
include(joinpath("gamma", "pions.jl"))


# datatypes and helper functions for LMB_SPECTRAL_CRs
export CRShockData,          # datatype to analyse single shocked particle
    readSingleCRShockDataFromOutputFile, # as the name says
    CRMomentumDistributionConfig, # config parameters for momentum distribution function
    CRMomentumDistribution,
    getCRMomentumDistributionFromPartID, # function to get distribution function
    get_detailled_shock_data,
    get_detailled_Dpp_data,
    get_detailled_radiative_data,
    get_detailled_adiabatic_data,
    synchrotron_emission,
    stokes_parameters,
    smallest_synch_bright_p,
    read_test_spectrum,
    construct_spectrum,
    momentum_bin_boundaries,
    momentum_bin_centers,
    get_synchrotron_spectrum_from_part_id,
    snapshot_to_spectra,
    write_cr_to_txt,
    read_cr_from_txt,
    write_cr_to_binary,
    read_cr_from_binary,
    cr_energy_in_range,
    energy_spectrum, 
    cr_number_in_range,
    gamma_source_PE04,
    gamma_source_pions,
    gamma_emissivity_pions,
    gamma_luminosity_pions,
    gamma_flux_pions


include("precompile.jl")
 
end # module
