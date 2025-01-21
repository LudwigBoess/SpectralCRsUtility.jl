```@meta
CurrentModule = SpectralCRsUtility
DocTestSetup = quote
    using SpectralCRsUtility
end
```

# Synchrotron emission

One important analysis for checking your simulations against observations is to compute the synchrotron emission of the simulated CR spectrum.
The can then be used to compute the total synchrotron power of the object, or to map the emission into a grid to construct mock observations.

You can compute the total synchrotron emission for a CR spectrum in a given magnetic field vie
```@docs
synchrotron_emission
```

If you are only interested in polarized emission you can compute this via
```@docs
synchrotron_polarisation
```

## Calculating Synchrotron emissivity

Please follow this example to compute the synchrotron emissivity at a given frequency for all SPH particles

```julia
using GadgetIO, GadgetUnits
using SpectralCRsUtility

filename = "path/to/file"

# read header
h = read_header(filename)

# define unit struct
GU = GadgetPhysical(h)

# read CR data
norm  = GU.CR_norm .* 10.0.^read_block(filename, "CReN", parttype=0)
slope = read_block(filename, "CReS", parttype=0) .|> Float64
cut   = read_block(filename, "CReC", parttype=0) .|> Float64

# cr setup 
pmin = 1.e-1           # minimum momentum defined in parameter file
pmax = 1.e5            # maximum momentum defined in parameter file
Nbins = size(slope,1)  # number of CR bins in the model

# define helper struct for synchrotron emission
par = CRMomentumDistributionConfig(pmin, pmax, Nbins)


# read magnetic field (if present)
bfld  = read_block(filename, "BFLD", parttype=0)
B = @. √( bfld[1,:]^2 + bfld[2,:]^2 + bfld[3,:]^2 )

# define observation frequency in Hz and shift to redshift
ν  = 1.44e9 * ( 1.0 + h.z )

# get total number of particles
Npart = length(cut)

# storage array for emissivities
j_ν = Vector{Float64}(undef, Npart)
for i = 1:Npart 
    j_ν[i] = synchrotron_emission(norm[:,i], slope[:,i], cut[i], B[i], par, ν0 = ν, 
                                    reduce_spectrum=true,
                                    integrate_pitch_angle = true)
end
```

# Gamma emission

Another form of emission from CRs is that of gamma-ray emission from pion-decay via inelastic scattering of CR protons with the thermal background protons.

This analysis is currently under construction and should be handled with care.

```@docs
gamma_source_pions
gamma_emissivity_pions
gamma_luminosity_pions
gamma_flux_pions
```