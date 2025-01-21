```@meta
CurrentModule = SpectralCRsUtility
DocTestSetup = quote
    using SpectralCRsUtility
end
```

# Basics

CRESCENDO sets up a distribution of piece-wise powerlaws to model a CR spectrum. 
The beginning and end of this powerlaw is defined in the OpenGadget3 parameter file. The number of bins are given in the Config file.

## Defining the spectrum

To define the CR spectrum for analysis with this package you can use
```@docs
CRMomentumDistributionConfig
```

## Momentum spectrum

To get an accurate representation of the simulated CR spectrum you need to reconstruct the spectrum from the Gadget output.
Gadget only stores the normalisation at the start of each momentum bin, the slope of the bin and the spectral cutoff.
In order to see if the solver is running fine you have to construct the normalisation at the end of each momentum bin from the slope of the bin.
If everything runs fine this should look similar to just plotting the normalisations. If there is a numerical issue it will show up as a stair-like feature.

You can construct the momentum distribution via
```@docs
CRMomentumDistribution
```
which returns a struct containing the bin boundaries and the normalisation at said boundaries.

## Number spectrum

To construct a spectrum from computing the number of CRs contained in each bin you can use
```@docs
energy_spectrum
```

# Energy contain in range

Another important step of analysis is to compute the total energy contained in the CR spectrum within a given energy band.

You can do this by using 
```@docs
cr_energy_in_range
```
