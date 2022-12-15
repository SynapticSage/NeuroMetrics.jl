"""
===========================
PFC FIRING DURING CA1 ISOLATION VERSUS ADJACENT
===========================
"""

using Infiltrator

function get_pos_labels(group; grouping, labels)
    p = []
    grouping = !(grouping isa Vector) ? [grouping] : grouping
    for g in grouping
        push!(p, Symbol(g) => group[1, g])
    end
    for (k,lab) in labels
        push!(p, Symbol(String(k) * "_label") => lab[group[1, k]])
    end
    #:cuemem => group.cuemem[1],
    #:cuemem_label => clab[group.cuemem[1]],
    p
end

function get_pfcrate_samples_at_spike(spikes, R; grouping=[], labels=Dict())
    pfc_rate_isoAdjacent_meanOfTimes = DataFrame()
    spikes = grouping == [] ? spikes : dropmissing(spikes, grouping)
    @time @showprogress for group in groupby(spikes, [:isolated,grouping...])
        #@info "samples" size(group.time)
        sample = [R[time=B, unit=At(pfc_units)]
                    for B in Between.(group.time.-0.015, group.time.+0.015)]
        sample = [s for s in sample if !isempty(s)]
        sample = vcat(sample...)
        unit = repeat(vec(pfc_units)', size(sample,1))
        append!(pfc_rate_isoAdjacent_meanOfTimes,
        DataFrame(OrderedDict(
            :unit => vec(unit),
            :rate => vec(sample),
            get_pos_labels(group; grouping=[:isolated, grouping...], labels)...
           )))
    end
    #@time pfc_rate_isoAdjacent_meanOfTimes.rankrate = sortperm(pfc_rate_isoAdjacent_meanOfTimes.rate)
    pfc_rate_isoAdjacent_meanOfTimes
end

function compute_differences_df(pfc_rate_isoAdjacent; 
        grouping=[], 
        addsamps=false,
        labels=Dict())
    D = DataFrame()
    grouping = Symbol.(grouping)
    groups = groupby(pfc_rate_isoAdjacent, grouping)
    @info "groups" length(groups)
    for group in groups
         adjacent_spikes = @subset(group,:isolated.==0).rate
         iso_spikes      = @subset(group,:isolated.==1).rate
         #radjacent_spikes = @subset(group,:isolated.==0).rankrate
         #riso_spikes      = @subset(group,:isolated.==1).rankrate
         zadjacent_spikes = zscore(@subset(group,:isolated.==0).rate)
         ziso_spikes      = zscore(@subset(group,:isolated.==1).rate)
         test = missing
         pval = missing
         rtest = missing
         rpval = missing
         try
             test =  UnequalVarianceTTest(iso_spikes, adjacent_spikes)
             pval = pvalue(test,tail=:right)
             #rtest =  UnequalVarianceTTest(riso_spikes, radjacent_spikes)
             #rpval = pvalue(rtest,tail=:right)
         catch
            continue
         end
         sampgroups = addsamps ? [
            :samp_iso => [iso_spikes],
            :samp_adj => [adjacent_spikes],
           ] : []
         
         append!(D, OrderedDict(
            :diff            => mean(iso_spikes) - mean(adjacent_spikes),
            #:rankdiff        => median(riso_spikes) - median(radjacent_spikes),
            :zdiff           => mean(ziso_spikes) - mean(zadjacent_spikes),
            :zdiv            => mean(ziso_spikes)/mean(zadjacent_spikes),
            :div             => mean(iso_spikes)/mean(adjacent_spikes),
            :iso_spikes      => mean(iso_spikes),
            :adjacent_spikes => mean(adjacent_spikes),
            :pval => pval,
            :test => test,
            :rpval => rpval,
            :rtest => rtest,
            get_pos_labels(group; grouping, labels)...,
            sampgroups...
           ))
    end
    D
end

"""
ARE PFC rates higher in isolated spiking?
"""

pfc_rate_isoAdjacent_meanOfTimes = get_pfcrate_samples_at_spike(spikes, R;
                                                                grouping=[],
                                                                labels=[])

#transform!(pfc_rate_isoAdjacent_meanOfTimes, DataFrames.All(), :isolated => (x->[isonames[xx] for xx in x]) => :isolated_label)
@subset!(pfc_rate_isoAdjacent_meanOfTimes, :isolated .== 0 .|| :isolated .== 1)
srt = UnequalVarianceTTest(@subset(pfc_rate_isoAdjacent_meanOfTimes,:isolated.==0).rate,
                           @subset(pfc_rate_isoAdjacent_meanOfTimes,:isolated.==1).rate)


# =========PLOTS =======================================
Plot.setfolder( "isolated_ca1_rate_pfc")
@df pfc_rate_isoAdjacent_meanOfTimes boxplot(:isolated, :rate, title="TTest=$(round(pvalue(srt),sigdigits=2)), with N=$(srt.df) PFC cells",
                         xticks=([0,1], collect(values(isonames))))
@df pfc_rate_isoAdjacent_meanOfTimes scatter!(:isolated, :rate, group=:isolated)
Plot.save((;save_kws...,desc="meanmean_pfc_cell_rate"))
# ======================================================

#"""
#===========================
#Mean of rates, per cell x time
#===========================
#"""
#
#D = compute_differences_df(pfc_rate_isoAdjacent_meanOfTimes; grouping=[], 
#                           labels=Dict(:isolated=>isonames))
#
#@subset!(D, :isolated .== 0 .|| :isolated .== 1)
#transform!(D, DataFrames.All(), :isolated => (x->[isonames[xx] for xx in x]) => :isolated_label)
#srt = UnequalVarianceTTest(@subset(D,:isolated.==0).rate,
#                          @subset(D,:isolated.==1).rate)


## =========PLOTS =======================================
#Plot.setfolder( "isolated_ca1_rate_pfc")
#@df isoadj_cellandtime boxplot(:isolated, :rate, title="Firing rate samples per spike\nTTest=$(round(pvalue(srt),sigdigits=2)) difference of individual cells",
#                         xticks=([0,1], collect(values(isonames))), outliers=false)
#Plot.save((;save_kws...,desc="firing rate samples per spike"))
## ======================================================
#
#"""
#===========================
#DO SAME split by cue and memory!
#===========================
#
#(STOPPED UPDATING HERE)
#"""
#
##TODO
#
#srt   = OrderedDict()
#diffs = OrderedDict()
#@showprogress for U in groupby(isoadj_cellandtime, :unit)
#    push!(srt,
#            U.unit[1] => UnequalVarianceTTest(@subset(U,:isolated.==0).rate,
#                                              @subset(U,:isolated.==1).rate))
#    push!(diffs,
#          U.unit[1] => mean(@subset(U,:isolated.==0).rate) - mean(@subset(U,:isolated.==1).rate))
#
#end
#
## =========PLOTS =======================================
#Plot.setfolder( "isolated_ca1_rate_pfc")
#pv_higheradjacent = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .> 0
#pv_higherisolated  = pvalue.(values(srt)) .< 0.05/length(srt) .&& values(diffs) .< 0
#pv_nonsig    = pvalue.(values(srt)) .> 0.05/length(srt)
#res =[collect(values(diffs)) pvalue.(values(srt))]
#bar([0, 1, 2], mean.([pv_nonsig, pv_higheradjacent, pv_higherisolated]),
#    xticks=([0, 1, 2], ["not bonferroni\nsig", "cells sig adjacent", "cells sig isolated"]),
#    ylabel="Percent", title="Nonlocality: PFC firing conditioned on CA1 field")
#
#Plot.save((;save_kws...,desc="bonferroni sig per cell, firing rate samples per spike"))
## ======================================================
#
#"""
#===========================
## Split by cuemem
#===========================
#"""
#
#Utils.filtreg.register(beh,spikes, on="time", transfer=["cuemem"])
#
#pfc_rate_isoAdjacent_meanOfTimes = get_pfcrate_samples_at_spike(spikes, R; 
#                            grouping=[:cuemem], labels=Dict(:cuemem => clab))
#D  = compute_differences_df(pfc_rate_isoAdjacent_meanOfTimes; 
#                            grouping=:cuemem, labels=Dict(:cuemem=>clab))
#D.clabel = getindex.([clab], D.cuemem)
#
#transform!(D, Colon(), :pval => Utils.pfunc => :pval_label)
##xticks=([-1,0,1], collect(values(clab)
#@df D bar(:clabel, :diff, group=:clabel)
#markers = text.(D.pval_label, :white)
#@df D annotate!(:clabel, :diff * 1.025, markers)
##transform!(Du, DataFrames.All(), :pval => Utils.pfunc => :pval_label)
#Plot.save((;save_kws..., desc="diff of iso adjacent in each state"))
#
#"""
#===========================
## Split by cuemem and pfc unit
#===========================
#"""
#
#pfc_rate_isoAdjacent_meanOfTimes = get_pfcrate_samples_at_spike(spikes, R; grouping=[:cuemem,:unit], labels=Dict(:cuemem=>clab))
#Du = compute_differences_df(pfc_rate_isoAdjacent_meanOfTimes; grouping=[:unit,:cuemem], labels=Dict(:cuemem=>clab))
#Du.sig = (Du.pval .< 0.05/21)
#Du.rsig = (Du.rpval .< 0.05/21)
#@df Du scatter(:cuemem + 0.1.*randn(size(:cuemem)), :diff, xticks=([-1,0,1], 
#    collect(values(clab))), group=:area, legend_position=:outerbottomright, 
#               ylabel="Δ(iso-adj) \$MUA_{pfc}\$ hz")
#hline!([0])
#
## SHEER NUMBER OF ISOLATED SPIKES (ALSO computer atop this script)
#if_counts = combine(groupby(Du, [:cuemem, :isolated]), nrow)
#subset!(if_counts, :isolated => x->(!).(isnan.(x)))
#if_perc = combine(groupby(if_counts, [:cuemem]), :isolated, :nrow => (x->x./(sum(x))) => :perc)
#iflab = Dict(0 => "out", 1=> "in")
#clab = Dict(-1 => "nontask", 0 => "cue", 1=> "mem")
#transform!(if_perc, DataFrames.All(), [:cuemem, :isolated] => 
#           ((x,y)->(getindex.([clab], x) .* " " .* getindex.([iflab], y))) .=> 
#           :cuemem_inoutfield)
#sort!(if_perc, :isolated)
#@df subset(if_perc, :isolated => (x->x.!=false)) bar(:cuemem_inoutfield, :perc, group=:isolated)
#
#Plot.save((;save_kws..., desc="cellcounts sig"))
#
#"""
#===========================
## PFC firign rates norm by 0-1
#===========================
#"""
#
#Plot.setfolder( "isolated_ca1_pfc_rates")
#Rpfc =  Munge.spiking.torate(@subset(allspikes, :area .== "PFC"), beh)
#Dpfc = DataFrame([0;Rpfc], :auto)
#Dpfc.cuemem = beh.cuemem
#pfc_cell_means = Matrix(combine(groupby(Dpfc, "cuemem"), 1:21 .=> mean))[:,2:end]
#pfc_cell_means./maximum(pfc_cell_means,dims=1)
#heatmap(["Nontask","Cue","Memory"], 1:size(pfc_cell_means,2), pfc_cell_means./maximum(pfc_cell_means,dims=1))
#Plot.save("Average PFC cell firing rates per task type")
#
##bar(["Nontask","Cue","Memory"], mean(pfc_cell_means./maximum(pfc_cell_means,dims=1), dims=2))
#
#Plot.save("Average of PFC cell firing, mean(norm0)")
#bar(["Nontask","Cue","Memory"], mean(pfc_cell_means, dims=2))
#Plot.save("Average of PFC cell firing, mean(rate)")
#
#"""
#===========================
## PFC firing rates norm by zscore
#===========================
#"""
#
#Plot.setfolder( "isolated_ca1_pfc_rates")
#Rpfc =  Munge.spiking.torate(@subset(allspikes, :area .== "PFC"), beh)
#Rpfc =  (Rpfc.-Matrix(mean(Rpfc,dims=1)))./Matrix(std(Rpfc,dims=1))
#Dpfc = DataFrame([0;Rpfc], :auto)
#Dpfc.cuemem = beh.cuemem
#pfc_cell_means = Matrix(combine(groupby(Dpfc, "cuemem"), 1:21 .=> mean))[:,2:end]
#pfc_cell_means./maximum(pfc_cell_means,dims=1)
#heatmap(["Nontask","Cue","Memory"], 1:size(pfc_cell_means,2), pfc_cell_means./maximum(pfc_cell_means,dims=1))
#Plot.save("Average PFC cell firing rates per task type")
#
#
#
#bar(["Nontask","Cue","Memory"], mean(pfc_cell_means, dims=2))
#Plot.save("Average of PFC cell firing, mean(norm0)")
#
#bar(["Nontask","Cue","Memory"], mean(pfc_cell_means, dims=2))
#Plot.save("Average of PFC cell firing, mean(rate)")

