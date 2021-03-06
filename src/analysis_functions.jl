
const cnst_c = 2.9979e10

@inline function energy_integral(bound_low::Real, bound_up::Real,
                                 norm::Real, slope::Real, ρ::Real)
    # energy integral (eq.21 M01)



    en = 4π * cnst_c * norm * bound_low^4 /  ρ
    energy = en * ( (bound_up/bound_low)^(4.0 - slope ) - 1.0 ) / ( 4.0 - slope )

    if ( 4.0 - slope_soft ) < slope < ( 4.0 + slope_soft )
        slope_var = (slope - 4.0)/slope_soft
        energy2 = en * log10(bound_up/bound_low)
        if slope_var != 0.0
            energy = energy * slope_var + energy2 * ( 1.0 - slope_var )
        else
            energy = energy2
        end
    end

    return energy
end

function calculateCREnergyInCGS(CR_N, CR_S, CR_Cut, ρ;
                                pmin::Real=10.0, pmax::Real=1.e7,
                                mc::Real=0.0, SelectBin::Int=-1,
                                units::GadgetPhysical=GadgetPhysical(),
                                verbose=true)
    # calculates engergy per bin in cgs units.

    if mc == 0.0
        error("Error! No mc specified! Please provide m_particle * c_light to continue!")
        return 1
    end

    CR_E = 0.0
    cnst_c = 2.9979e10
    ρ *= units.m_unit / units.l_unit^3  # convert rho in cgs units

    CR_N .*= 1.e20

    Nbins = length(CR_N)
    bin_width = log10(pmax/pmin)/Nbins
    bin = 10.0.^collect(log10(pmin):bin_width:log10(pmax)) .* mc
    above_cut = findall( bin .> CR_Cut/mc )
    if length(above_cut) > 0
        bin[above_cut] .= CR_Cut/mc
    end

    if SelectBin == -1
        CR_E = zeros(length(CR_N))
        for i = 1:length(CR_N)
            CR_E[i] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    elseif length(SelectBin) == 1
        CR_E = energy_integral(bin[SelectBin], bin[SelectBin+1], CR_N[SelectBin], CR_S[SelectBin], ρ)
    elseif length(SelectBin) > 1
        CR_E = zeros(length(SelectBin))
        if i ∈ SelectBin
            CR_E[i-SelectBin[1]+1] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    end
    return CR_E
end


@inline function density_integral(bound_low::Real, bound_up::Real,
                         norm::Real, slope::Real, ρ::Real)
    # density integral (eq. 9 M01)

    nb = 4π * norm * bound_low^3 / ρ
    density = nb * ( (bound_up/bound_low)^(3.0 - slope ) - 1.0 ) / ( 3.0 - slope )

    if ( 3.0 - slope_soft ) < slope < ( 3.0 + slope_soft )
        slope_var = (slope - 3.0)/slope_soft
        density2 = nb * log10(bound_up/bound_low)
        if slope_var != 0.0
            density = density * slope_var + density2 * ( 1.0 - slope_var )
        else
            density = density2
        end
    end

    return density
end


function calculateCRNumber(CR_N, CR_S, CR_Cut, ρ;pmin::Real=10.0, pmax::Real=1.e7,
                           mc::Real=0.0, SelectBin::Int=-1,
                           units::GadgetPhysical=GadgetPhysical(), verbose=true)
    # calculates engergy per bin in cgs units.

    CR_Num = 0.0
    cnst_c = 2.9979e10
    ρ *= units.m_unit / units.l_unit^3  # convert rho in cgs units

    CR_N .*= 1.e20  # compensate for io

    Nbins = length(CR_N)
    bin_width = log10(pmax/pmin)/Nbins
    bins = 10.0.^collect(log10(pmin):bin_width:log10(pmax)) .* mc
    above_cut = findall( bins .> CR_Cut/mc )

    if length(above_cut) > 0
        bins[above_cut] .= CR_Cut/mc
    end

    if SelectBin == -1
        CR_Num = zeros(length(CR_N))
        for i = 1:length(CR_N)
            CR_Num[i] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    elseif length(SelectBin) == 1
        CR_Num = energy_integral(bin[SelectBin], bin[SelectBin+1], CR_N[SelectBin], CR_S[SelectBin], ρ)
    elseif length(SelectBin) > 1
        CR_Num = zeros(length(SelectBin))
        if i ∈ SelectBin
            CR_Num[i-SelectBin[1]+1] = energy_integral(bin[i], bin[i+1], CR_N[i], CR_S[i], ρ)
        end
    end
    return CR_Num
end
