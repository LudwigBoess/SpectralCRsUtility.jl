import PyPlot.plot 

plot(cr::CRMomentumDistribution; kwrdargs...) = plot(cr.bound[1:end-1], cr.norm; kwrdargs...)
