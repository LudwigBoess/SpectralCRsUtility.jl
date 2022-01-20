[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md)
[![Build Status](https://travis-ci.org/LudwigBoess/SpectralCRsUtility.jl.svg?branch=master)](https://travis-ci.org/LudwigBoess/SpectralCRsUtility.jl)
[![codecov.io](https://codecov.io/gh/LudwigBoess/SpectralCRsUtility.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/SpectralCRsUtility.jl?branch=master)

# SpectralCRsUtility.jl

Provides helper functions and analysis function for LMB_SPECTRAL_CRs in OpenGadget3 (not public).

## Calculating Synchrotron emissivity

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