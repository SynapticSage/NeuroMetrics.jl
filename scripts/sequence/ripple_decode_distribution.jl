import Plots

ca1, ca1pfc, pfcca1, pfc = ripple[:,:,:,1], 
                           ripple[:,:,:,2], 
                           ripple[:,:,:,3],
                           ripple[:,:,:,4]
nca1, nca1pfc, npfcca1, npfc = (!).(isnan(ripple[:,:,:,1])), 
   (!).(isnan(ripple[:,:,:,2], 
   (!).(isnan(ripple[:,:,:,3],
   (!).(isnan(ripple[:,:,:,4]

mean_dist = Dict()
mean_dist[:pfc] = utils.squeeze(nanmean(pfc; dims=(3)))
mean_dist[:ca1] = utils.squeeze(nanmean(ca1; dims=(3)))
mean_dist[:pfcca1] = utils.squeeze(nanmean(pfcca1; dims=(3)))
mean_dist[:ca1pfc] = utils.squeeze(nanmean(ca1pfc; dims=(3)))

Plots.plot(
Plots.heatmap(mean_dist[:ca1],title="CA1"),
Plots.heatmap(mean_dist[:pfc],title="PFC"),
Plots.heatmap(mean_dist[:ca1pfc],title="CA1PFC"),
Plots.heatmap(mean_dist[:pfcca1],title="PFCCA1")
)
Plots.plot(
Plots.plot(Plots.heatmap(mean_dist[:pfc]-mean_dist[:ca1], c=:vik, title="PFC > CA1")),
Plots.heatmap(mean_dist[:pfcca1]-mean_dist[:ca1pfc], c=:vik, title="PFCCA1 > CA1PFC"),
Plots.heatmap(mean_dist[:ca1pfc]-mean_dist[:pfcca1], c=:vik, title="CA1PFC > PFCCA1")
)
