"""
        INFO HERE
"""

"""
                Debug data from shock
"""
mutable struct CRShockData
   """   Data structure to analyze a single shocked particle.
         To be used with compile option: LMB_CR_DEBUG_TRACK_PARTICLE
   """
   dt::Vector{Float64}
   Mach::Vector{Float64}
   Shock_Speed::Vector{Float64}
   Shock_Compress::Vector{Float64}
   Shock_Energy_In::Vector{Float64}
   Shock_Energy_Real::Vector{Float64}
   Energy_P::Vector{Float64}
   Energy_e::Vector{Float64}

   function CRShockData(N::Int64=0)

       new(zeros(N), zeros(N), zeros(N), zeros(N),
           zeros(N), zeros(N), zeros(N), zeros(N))

   end
end



"""
                Distribution spectrum in momentum space
"""

