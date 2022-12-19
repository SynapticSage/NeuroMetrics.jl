quickactivate(expanduser("~/Projects/goal-code/"));
using GoalFetchAnalysis
using Timeshift
using Timeshift.types
using Timeshift.shiftmetrics
using Field.metrics
using Plot
using Plot.receptivefield
using Utils.namedtup
using Munge.timeshift: getshift
using Munge.nonlocal
using Utils.statistic: pfunc
import Plot
using Munge.spiking
using Filt
using Infiltrator
using SoftGlobalScope

using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
using LazySets
GC.gc()
datasets = (("RY16",36, nothing),("RY22",21, ""), ("RY22",21, "ref"))
for (animal,day,ref) in datasets

    clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
    Munge.nonlocal.setclab(clab)
    isonames =  OrderedDict(false => :adjacent, true=>:isolated)
    filt_desc = OrderedDict(:all => "> 2cm/s")
    save_kws = (;pfc_rate_analy=true)
    filt = Filt.get_filters()
    datacut = :all

    Plot.setappend((;animal,day,ref))
    Plot.setparentfolder("nonlocality")

    # ===================
    # ACQUIRE DATA
    # ===================
    # Acquire data
    @time spikes, beh, cells = Load.load(animal, day, data_source=["spikes","behavior", "cells"])
    GC.gc()
    beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], 
    filters=filt[datacut], filter_skipmissingcols=true)
    allspikes = copy(spikes)
    beh2 = Load.load_behavior(animal,day)
    Munge.nonlocal.setunfilteredbeh(beh2)

    # THis section requires enormous RAM, 100GB at least
    LFP = Load.load_lfp(animal, day; ref)
    ca1 = sort(unique(@subset(cells, :area .== "CA1").tetrode))
    LFP = groupby(LFP, :tetrode)
    LFP = [LFP[(;tetrode)] for tetrode in ca1]
    GC.gc()
    LFP = vcat(LFP...)
    GC.gc()

    using SoftGlobalScope

    iso, phase = [], []
    @softscope for tet in ca1

        @info "loop start" tet length(iso) length(phase)

        lfp = copy(@subset(LFP, :tetrode .== tet))
        lfp.time = lfp.time .- Load.min_time_records[end]
        try
            lfp = Munge.lfp.annotate_cycles(lfp, method="peak-to-peak") # TODO potential bug, 1st time runs, cuts trough-to-trough, second peak-to-peak
        catch
            push!(iso,nothing)
            push!(phase,nothing)
            continue
        end
        @assert length(unique(lfp.cycle)) > 1
        #sp = @subset(spikes, :tetrode .== 6);
        @df lfp[1:2500,:] begin
            Plots.plot(:time, :raw, label="raw")
            Plots.plot!(:time, mod2pi.(:phase) .+100, label="phase")
            Plots.plot!(:time, 10*:cycle)
        end
        Utils.filtreg.register(lfp, spikes, on="time", transfer=["phase"])
        @assert !all(ismissing.(spikes.phase))

        # ===================
        # ISOLATED SPIKING
        # ===================
        #Munge.spiking.isolated(sp, lfp, include_samples=true)
        Munge.spiking.isolated(spikes, lfp, include_samples=false, 
                               refreshcyc=true, overwrite=true)
        push!(iso, spikes.isolated)
        push!(phase, spikes.phase)
    end

    import JLD2
    JLD2.save(datadir("checkpoint_tetrodecheck_$(animal)_$(day).tmp.jld2"), "iso", "phase")

    df = vcat((DataFrame(Dict("tetrode"=>tet,"phase"=>ph,"isolated"=>is)) for (tet,ph,is) in zip(ca1, phase, iso) 
               if ph !== nothing)...)
    gdf = groupby(df, :tetrode)
    pl = [(@df DataFrame(dropmissing(g, :isolated)) histogram(:phase, group=:isolated, bins=40; normalize=:probability,label="",xticks=[],yticks=[], alpha=0.4, title="$(g.tetrode[1])")) for g in gdf]
    plot(pl...)

    Plot.setparentfolder("isolated","find_goot_tet.jl")
    Plot.save("tetrodePhaseLocking")

end
