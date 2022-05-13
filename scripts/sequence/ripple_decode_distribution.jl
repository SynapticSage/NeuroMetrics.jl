quickactivate("/home/ryoung/Projects/goal-code")
include(scriptsdir("decode","Initialize.jl"))
include(scriptsdir("decode","LoadData.jl"))
import Plots

# TODO list
# Retry this with unfiltered dat, no quantile filtering, and generate ripple
# array structure

D = decode.load_checkpoints(decode_file, vars=[:ripple])

ca1, ca1pfc, pfcca1, pfc = view(ripple,:,:,:,1), 
                           view(ripple,:,:,:,2), 
                           view(ripple,:,:,:,3),
                           view(ripple,:,:,:,4)

nca1, nca1pfc, npfcca1, npfc = (!).(isnan.(ripple[:,:,:,1])), 
                                (!).(isnan.(ripple[:,:,:,2])), 
                                (!).(isnan.(ripple[:,:,:,3])),
                                (!).(isnan.(ripple[:,:,:,4]))

mean_dist = Dict()
mean_dist[:pfc] = utils.squeeze(nanmean(pfc; dims=(3)));
mean_dist[:ca1] = utils.squeeze(nanmean(ca1; dims=(3)));
mean_dist[:pfcca1] = utils.squeeze(nanmean(pfcca1; dims=(3)));
mean_dist[:ca1pfc] = utils.squeeze(nanmean(ca1pfc; dims=(3)));
nmean_dist = Dict()
nmean_dist[:pfc] = utils.squeeze(nanmean(npfc; dims=(3)));
nmean_dist[:ca1] = utils.squeeze(nanmean(nca1; dims=(3)));
nmean_dist[:pfcca1] = utils.squeeze(nanmean(npfcca1; dims=(3)));
nmean_dist[:ca1pfc] = utils.squeeze(nanmean(nca1pfc; dims=(3)));
nmedian_dist = Dict()
nmedian_dist[:pfc]    = utils.squeeze(nanmedian(npfc; dims=(3)));
nmedian_dist[:ca1]    = utils.squeeze(nanmedian(nca1; dims=(3)));
nmedian_dist[:pfcca1] = utils.squeeze(nanmedian(npfcca1; dims=(3)));
nmedian_dist[:ca1pfc] = utils.squeeze(nanmedian(nca1pfc; dims=(3)));



Plots.plot(
Plots.heatmap(mean_dist[:ca1],title="CA1", clim=(0, 0.07)),
Plots.heatmap(mean_dist[:pfc],title="PFC", clim=(0, 0.07)),
Plots.heatmap(mean_dist[:ca1pfc],title="CA1PFC", clim=(0, 0.07)),
Plots.heatmap(mean_dist[:pfcca1],title="PFCCA1", clim=(0, 0.07)),
)
Plots.plot(
Plots.plot(Plots.heatmap(mean_dist[:pfc]-mean_dist[:ca1], c=:vik, title="PFC > CA1")),
Plots.heatmap(mean_dist[:pfcca1]-mean_dist[:ca1pfc],      c=:vik, title="PFCCA1 > CA1PFC"),
Plots.heatmap(mean_dist[:ca1pfc]-mean_dist[:pfcca1],      c=:vik, title="CA1PFC > PFCCA1")
)

p = Plots.plot(
Plots.heatmap(nmean_dist[:ca1],title="CA1"),
Plots.heatmap(nmean_dist[:pfc],title="PFC"),
Plots.heatmap(nmean_dist[:ca1pfc],title="CA1PFC"),
Plots.heatmap(nmean_dist[:pfcca1],title="PFCCA1"),
size=(1200,800),
)
utils.savef(plotsdir("ripples", "ca1,pfc,coordinated", "mean_counts_above_97quantile"))

Plots.plot(
Plots.heatmap(nmedian_dist[:ca1],    title="CA1"),
Plots.heatmap(nmedian_dist[:pfc],    title="PFC"),
Plots.heatmap(nmedian_dist[:ca1pfc], title="CA1PFC"),
Plots.heatmap(nmedian_dist[:pfcca1], title="PFCCA1"),
)
