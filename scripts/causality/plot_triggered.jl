#                    
#,---.|         |    
#|---'|    ,---.|--- 
#|    |    |   ||    
#`    `---'`---'`---'
#
# Plots the results of a triggered causal run
# where the run has x-y variables in its split
#
# Before running this, you need data in your workspace
using DrWatson
using Plot.task
include(scriptsdir("causality", "init_trig_plot.jl"))

func = nanmean
tsk = Load.load_task(animal, day)

function getnanfrac(X)
    horizonAndSamp = Tuple(collect(1:ndims(X))[end-1:end])
    total_who_could_be_present = size(X,horizonAndSamp[1]) .* 
                                 size(X,horizonAndSamp[2])
    sum(isnan.(X); dims=horizonAndSamp) ./
                        total_who_could_be_present
end
"""
    plottrackcause(ca1pfc,pfcca1, 
                     slice=[Colon() for _ in 1:ndims(ca1pfc)], 
                          plotdims=1:2; weightstyle=:nanfrac)

Plot the causal data in ca1pfc and pfcca1 matrices, indexed by `slice`, and
along the plotdimensions `plotdims`
"""
function plottrackcause(ca1pfc,pfcca1, 
                          slice=[Colon() for _ in 1:ndims(ca1pfc)], 
                          plotdims=1:2; weightstyle=:nanfrac)

    axes = [collect(Float64.(c)) for c in grd.centers][plotdims]

    # How many unsampled?
    nanfrac = getnanfrac(ca1pfc)

    # Perform slicing
    ca1pfc = ca1pfc[slice...]
    pfcca1 = pfcca1[slice...]
    nanfrac = nanfrac[slice...]

    squishdims = (setdiff(1:ndims(ca1pfc), collect(plotdims))...,)
    ca1pfc  = nanmean(ca1pfc,dims=squishdims)
    pfcca1  = nanmean(pfcca1;dims=squishdims)
    nanfrac = nansum(nanfrac;dims=squishdims)

    # Make sure the final product is XD
    ca1pfc  = Utils.squeeze(ca1pfc)
    pfcca1  = Utils.squeeze(pfcca1)
    nanfrac = Utils.squeeze(nanfrac)

    if weightstyle == :nanfrac
        weighting = (1 .- nanfrac)
    elseif weightstyle == :sampcount
        weighting = weighting_trig'
    elseif weightstyle == :nanfrac_sq
        weighting = ((1 .- nanfrac) .^ 2)
    end
    weigthing = replace(weighting, 0.0 => NaN)

    pfcca1  = (pfcca1 .* weighting)[:,:]'
    ca1pfc  = (ca1pfc .* weighting)[:,:]'
    nanfrac = nanfrac'

    @assert ndims(ca1pfc)<=2 "Dafuq. Dims need to be no more than 2."
    
    p_unsamp=heatmap(axes..., nanfrac, title="unsampled fraction")
    plotboundary!(tsk,transpose=true, c=:black)

    
    h_pfcca1 = heatmap(axes[1], axes[2], pfcca1, title="pfc-ca1", 
                     c=:vik, clims=nanmaximum(abs.(pfcca1)).*(-1,1))
    plotboundary!(tsk,transpose=true, c=:black)

    h_ca1pfc=heatmap(axes[1], axes[2], ca1pfc, title="ca1-pfc", 
                     c=:vik, clims=nanmaximum(abs.(ca1pfc)).*(-1,1))
    plotboundary!(tsk,transpose=true, c=:black)

    plot(h_pfcca1,h_ca1pfc, size=(800,400))
    plot(h_pfcca1,h_ca1pfc,p_unsamp, size=(800,800))

    dif = pfcca1 .- ca1pfc
    dif = func(dif .* weighting, dims=(3,4))[:,:]' .* (2)
    hD = heatmap(axes..., dif; title="pfcca1 .- ca1pfc", c=:vik, clims=nanmaximum(abs.(dif)).*(-1,1))
    plotboundary!(tsk,transpose=true,c=:black,linestyle=:solid,label="")

    plot(h_pfcca1,h_ca1pfc,p_unsamp,hD, size=(800,400))
end
@time plottrackcause()

anim = @animate for i in params[:horizon]
    plottrackcause(i)
end
gif(anim, fps=5)
