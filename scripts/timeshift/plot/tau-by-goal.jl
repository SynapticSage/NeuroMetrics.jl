using Serialization
using ProgressMeter
import Plot
using DataStructures: OrderedDict
using Plots
using HypothesisTests
using SoftGlobalScope
using DataFrames
using Timeshift
using Timeshift.shiftmetrics
import Field.metrics
using Statistics, NaNStatistics, HypothesisTests
using DimensionalData
filts = Filt.get_filters_precache()
thresh = 3f0
shifts = -2f0:0.05f0:2f0

#animal,day="RY22",21
@softscope for (animal,day) in (("super",0),) # ("RY16", 36))
    @info "dataset" animal day
    @time spikes, beh, ripples, cells = Load.load(animal, day)
    Plot.setappend("$animal-$day")

    w = 8
    Plot.setfolder("timeshift","xyG-xywidth=$w")

    props = ["x","y","stopWell"]
    widths = OrderedDict(k=>v for (k,v) in zip(props, Field.getwidths(props)))
    widths["stopWell"] = 0.50f0
    widths["x"] = w
    widths["y"] = w
    widths = collect(values(widths))
    Utils.filtreg.register(beh,spikes,on="time",transfer=["stopWell"])

    # NOTE :
    # IF you want to use :nonatsk in here, you need to have pseudo-stopWell
    # determined. This means using state-transition history to determine
    # the best fiting route @ a given time

    #@showprogress for datacut in [:all, :task, :cue, :memory]
    @showprogress for datacut in [:all, :task, :cue, :memory, :correct, :error]
        @info datacut
        global shifted = Timeshift.shifted_fields(dropmissing(beh,:stopWell),
                                       dropmissing(spikes, :stopWell), 
                                       shifts, ["x","y"];
                                       shiftbeh=false,
                                       widths, 
                                       adaptive=false,
                                       metricfuncs=[metrics.bitsperspike,
                                                    metrics.totalcount,
                                                    metrics.maxrate,
                                                    metrics.meanrate],
                                       filters=filts[datacut], 
                                       thresh);

        global shifted_wg = Timeshift.shifted_fields(dropmissing(beh,:stopWell),
                                       dropmissing(spikes, :stopWell), 
                                       shifts, ["x","y","stopWell"];
                                       shiftbeh=false,
                                       widths, 
                                       adaptive=false,
                                       metricfuncs=[metrics.bitsperspike,metrics.totalcount,metrics.maxrate,metrics.meanrate],
                                       filters=filts[datacut], 
                                       thresh);

        fwg, f = matrixform(ShiftedFields(shifted_wg)), matrixform(ShiftedFields(shifted))

        serialize(datadir("exp_pro", "$animal-$day-xyG-$datacut-$w-fmat"), (;f, fwg, animal, day, shifts, thresh, w))

        # Get shifts
        sh = collect(shifts)

        # Add metrics
        push_shiftmetric!(fwg, best_tau!; metric=:bitsperspike)

        # -------------
        # GOAL INDEXING
        # -------------
        # Grab the goalindex as the mean response to goal-pursuit across the
        # different XY positions 
        gd = [Utils.squeeze(nanmean(nanmean(f.rate; dims=1), dims=2))
              for f in fwg]
        gm = [nanmaximum(g)-nanminimum(g) for g in gd]
        g_count = [Utils.squeeze(nansum(nansum(f.count; dims=1), dims=2))
                    for f in fwg]
        g_occ_count = [Utils.squeeze(nansum(nansum(f.occ.count; dims=1), dims=2))
                    for f in fwg]
        gc = [nanmaximum(c./o)-nanminimum(c./o)
              for (c,o) in zip(g_count, g_occ_count)]

        # CLEAN :  Throw away responses with best BPS less than 0.5 (papers throw out these neurons)
        fwg[:bestshift_bitsperspike][fwg[:bestshift_bitsperspike] .< 0.5] .= NaN

        # Obtain a dimarray representation
        gm = DimArray(gm, fwg.dims)

        B = []
        for (r_fwg, r_gm) in zip(eachrow(fwg),eachrow(gm))
            samp = [r_fwg[shift=At(r_fwg[:bestshift_bitsperspike][1])][:bestshift_bitsperspike],
             r_gm[shift=At(r_fwg[:bestshift_bitsperspike][1])]]
            push!(B, samp)
        end
        B = hcat(B...)'
        B[:,1] = abs.(B[:,1])


        se(x)  = std(x)/sqrt(length(x))

        inds = B[:,1] .> 0.25
        chunks = [(B[(!).(inds),2]), (B[inds,2])]
        t = OneWayANOVATest(chunks...)
        bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")
        Plot.save("summary <0,25, datacut=$datacut, w=$w, pval=$(pvalue(t))")
        
        inds = B[:,1] .> 0.5
        chunks = [(B[(!).(inds),2]), (B[inds,2])]
        t = OneWayANOVATest(chunks...)
        bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")
        Plot.save("summary, datacut=$datacut, w=$w, pval=$(pvalue(t))")

        inds = B[:,1] .> 1
        chunks = [(B[(!).(inds),2]), (B[inds,2])]
        t = OneWayANOVATest(chunks...)
        bar(1:2, mean.(chunks), yerror=se.(chunks), title="pval=$(pvalue(t))")
        Plot.save("summary, <1, datacut=$datacut, w=$w, pval=$(pvalue(t))")
        ##scatter(vec(fwg[:bestshift_bitsperspike]), vec(gm), title="$datacut")
 
        scatter(vec(B[:,1]), vec(B[:,2]), title="$datacut")
        Plot.save((;desc="scatter goal_index versus bestshift_bitsperspike",
                   datacut))

        histogram2d(vec(fwg[:bestshift_bitsperspike]), vec(gm), title="$datacut")
        Plot.save((;desc="histogram2d goal_index versus bestshift_bitsperspike",
                   datacut))

        histogram(vec(fwg[:bestshift_bitsperspike]), title="$datacut")
        Plot.save((;desc="histogram bestshift_bitsperspike", datacut))
        inds = sortperm(fwg[:bestshift_bitsperspike][:,1])
        bps  = fwg[:bitsperspike][inds, :]

        #X = hcat([Utils.norm_extrema(b) for b in eachrow(bps)]...)'
        #heatmap(sh, collect(1:size(X,1)), X, title="$datacut", size=(400,1000))
        #vline!([0],c=:black,linestyle=:dash, linewidth=2)
        #Plot.save((;desc="snake plot, norm 01", datacut))

        #X = hcat([Utils.norm_percent(b,0.5) for b in eachrow(bps)]...)'
        #heatmap(sh, collect(1:size(X,1)), X, title="$datacut", size=(400,1000))
        #vline!([0],c=:black,linestyle=:dash, linewidth=2)
        #Plot.save((;desc="snake plot, norm percent", datacut))
    end
end

#fwg[:sortshift_bps] # TODO create this
#scatter(vec(gm), vec(fwg[:sortshift_bps]))
