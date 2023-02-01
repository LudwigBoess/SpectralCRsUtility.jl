using QuadGK


const j_nu_prefac_miniati = √(3)*q_e^3 / (2cL)
const νB_fac = q_e / ( 2π*m_e*cL )

F(x)  = 0.87*∛(x) * exp( -11/8 * x^(7/8) )
get_νB(B) = νB_fac * B

𝒳(p, ν, νB) = 1.5 * (ν/νB) / p^2


function synchrotron_emission_miniati(  f_p::Vector{<:Real},
                                        q::Vector{<:Real},
                                        cut::Real,
                                        B_cgs::Real,
                                        bounds::Vector{<:Real};
                                        ν0::Real = 1.4e9,
                                        integrate_pitch_angle::Bool = true,
                                        convert_to_mJy::Bool = false,
                                        reduce_spectrum::Bool = true)


    j_ν = 0.0

    νB = get_νB(B_cgs)

    for i = 1:Nbins

        # if the lower boundary is above the cutoff all other bins don't contribute
        if bounds[i] > cut 
            break
        end

        # default upper boundary
        p_up = bounds[i+1]

        # upper boundary if cutoff is within the bin
        if p_up > cut 
            p_up = cut
        end

        # helper function 
        f(p) = F(𝒳(p, ν0, νB)) * 𝒳(p, ν0, νB)^(-0.5 * (q[i] - 5))

        # x-integral
        integral, err = quadgk(f, bounds[i], p_up, rtol=1e-8)

        # sum up bin contribution
        j_ν += f_p[i] * (bounds[i]*m_e*cL)^q[i] * (2ν0 / 3νB)^( -0.5 * (q[i] - 3) ) * integral
    end

    return j_ν * j_nu_prefac_miniati
end


using GadgetUnits
using BenchmarkTools
using SpectralCRsUtility

GU = GadgetPhysical()
GU.CR_norm

pmin = 1.0
pmax = 1.e6
Nbins = length(ref_cr_norm)
bin_width = log10(pmax/pmin)/Nbins
bounds = [pmin * 10.0^((i - 1) * bin_width) for i = 1:Nbins+1]

jν_miniati = synchrotron_emission_miniati(GU.CR_norm .* ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, bounds)

jν = synchrotron_emission(GU.CR_norm .* ref_cr_norm, ref_cr_slope, ref_cr_cut, ref_cr_B, bounds)

norm = GU.CR_norm .* ref_cr_norm
@btime synchrotron_emission_miniati($norm, $ref_cr_slope, $ref_cr_cut, $ref_cr_B, $bounds)

@btime synchrotron_emission($norm, $ref_cr_slope, $ref_cr_cut, $ref_cr_B, $bounds)
