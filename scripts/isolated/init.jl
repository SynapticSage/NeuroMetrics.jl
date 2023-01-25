using DrWatson
quickactivate(expanduser("~/Projects/goal-code/"));

using GoalFetchAnalysis
using Timeshift, Timeshift.types, Timeshift.shiftmetrics
using Field.metrics
using Plot, Plot.receptivefield
using Munge.spiking, Munge.nonlocal
using Munge.timeshift: getshift
using Utils.namedtup 
using Utils.statistic: pfunc
using Filt

using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
using LazySets
using JLD2
using Infiltrator

datasets = (("RY16",36,:ca1ref), ("RY22",21,:ca1ref),  #("super", 0, :ca1ref),
            ("RY16",36,:default),("RY22",21,:default), #("super", 0, :default)
           )
(animal, day, tet) = datasets[2]

opt = Dict(:process_outoffield => true, :ploton=>true)

# Plotting turned on?
if opt[:ploton] 
    Plot.on()
else
    Plot.off()
end

init = 1
if init != 1; @warn("initial dataset is $init"); end
(animal, day, tet) = datasets[1]
@showprogress "datasets" for (animal, day, tet) in datasets[init:end]

    @info "loop" animal day tet

    clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
    Munge.nonlocal.setclab(clab)
    isonames =  OrderedDict(false => :adjacent, true=>:isolated)
    filt_desc = OrderedDict(:all => "> 2cm/s")
    save_kws = (;pfc_rate_analy=true)
    filt = Filt.get_filters()
    datacut = :all

    Plot.setappend((;animal,day,tet))

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
    lfp = Load.load_lfp(animal, day, tet=tet);

    # =============
    ## MUNGE LFP ##
    # =============
    # Center and annotate
    begin
        lfp.time = lfp.time .- Load.min_time_records[end]
        # TODO potential bug, 1st time runs, cuts trough-to-trough, second peak-to-peak
        lfp = Munge.lfp.annotate_cycles(lfp, method="peak-to-peak") 
    end
    # Visualize our annotations
    begin
        @assert length(unique(lfp.cycle)) > 1
        #sp = @subset(spikes, :tetrode .== 6);
        @df lfp[1:2500,:] begin
            Plots.plot(:time, :raw, label="raw")
            Plots.plot!(:time, mod2pi.(:phase) .+100,label="phase")
            Plots.plot!(:time, 10*:cycle, label="cycle labels")
        end
    end

    # Transfer lfp phase to spikes (for phase locking measurements)
    begin
        Utils.filtreg.register(lfp, spikes, on="time", transfer=["phase"])
        @assert !all(ismissing.(spikes.phase))
    end
    ##


    # ===================
    # ISOLATED SPIKING
    # ===================
    #Munge.spiking.isolated(sp, lfp, include_samples=true)
    Munge.spiking.isolated(spikes, lfp, include_samples=false)

    # ===================
    # FIELDS
    # ===================
    if opt[:process_outoffield]
        F = load_fields()
        kz = collect(filter(k->k.animal == animal && k.day == day, keys(F)))
        @time f = F[bestpartialmatch(kz, 
                                     (;datacut, widths=5,coactivity=nothing), 
                                     nothing_means_removekey=true)];
        f = f isa ShiftedFields ? f : ShiftedFields(deepcopy(f))
        unitshift = Timeshift.types.matrixform(f)

        # ===================
        # OUT OF FIELD SPIKES
        # ===================
        # Setup a  shift-getting convenience method, the shifts, and a few metrics
        shifts = collect(unitshift.dims[2])
        push_dims!(unitshift)
        push_celltable!( unitshift, cells, :unit, :area)
        push_metric!(unitshift, Field.metrics.bitsperspike)
        push_metric!(unitshift, Field.metrics.totalcount)
        push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)
        annotate_nonlocal_spikes!(spikes, cells, unitshift, 0)
        # Which cells pass our criteria?
        #region = :CA1
        #metricfilter = metricfilters[region]
        Plot.setparentfolder("nonlocality")
    end


    @assert length(unique(lfp.cycle)) > 1
    cycles = Munge.lfp.get_cycle_table(lfp)
    Munge.lfp.annotate_cycles!(spikes, cycles)

    # ENSURING NO SLEEP :: Usually taken care of
    begin
        tsk = Load.load_task(animal, day)
        tsk = DataFrame(unique(eachrow(tsk[:,[:start,:end,:task]])))
        tsk = subset(tsk, :task=>t->t .!= "sleep")
        ismissing.(spikes.velVec)
    end

    # SAVE THE DATA
    begin
        filename = datadir("isolated","iso_animal=$(animal)_day=$(day)_tet=$(tet).jld2")
        !isdir(dirname(filename)) ? mkpath(dirname(filename)) : nothing
        @save "$filename" {compress=true} animal day lfp spikes allspikes tsk cells beh cycles
    end

    Plot.setparentfolder("isolated")

    # Testing isolation
    Plot.setfolder("phase_locking")
    # All spikes, no cell grouping
    begin
        spikes = :phase ∈ propertynames(spikes) ? spikes[!,Not(:phase)] : spikes
        Utils.filtreg.register(lfp, spikes, on="time", transfer=["phase"])
        SP = subset(dropmissing(spikes,:isolated),:velVec => v->abs.(v) .> 3)
        SP.pyrint = replace(SP.meanrate .> 4, 0=>"pyr",1=>"int")
        P=[]
        for sp in groupby(SP, :area)
            s = subset(sp, :pyrint => s->s.=="pyr")
            p1=histogram(s.phase, group=s.isolated, normalize=:pdf, alpha=0.5, 
                      legend_title=:isolated, ylabel="fraction",
                      title="PL of $(sp.area[1]) to $tet tetrode"
                     )
            lfp.phase_digitize = Int8.(Utils.binning.digitize(lfp.phase, 80+1))
            p2=@df combine(groupby(lfp, :phase_digitize), 
                           :phase=>x->maximum(abs.(x)).*mode(sign.(x)),
                    :raw => mean, renamecols=false) plot(:phase, :raw; xlabel="phase", c=:black, linewidth=2)
            lay = @layout [a{0.8h}; b{0.2h}]
            p=plot(p1,p2, layout=lay, label="")
            push!(P, p)
            p
            Plot.setfolder("phase_locking")
            Plot.save("all cell in $(sp.area[1]) to $tet tetrode")
        end
    end
    P

    Plot.setfolder("phase_locking", "cells_$animal-$day-$tet")
    begin # Individual cells
        sp = subset(dropmissing(spikes,:isolated),
                        :velVec => v->abs.(v) .> 4)
        sp.pyrint = replace(sp.meanrate .> 5, 0=>"pyr",1=>"int")
        for sp_area in groupby(sp, :area)

            @info sp_area.area[1]
            H = []
            G=  groupby(sp_area, :unit)
            @showprogress for s in G
                title = if s === first(G)
                    "$area=$(s.area[1])" 
                else
                    ""
                end
                h=histogram(s.phase; 
                            group=s.isolated, normalize=:pdf, bins=20, title,
                            alpha=0.5, legend_title=:isolated, xlabel="phase",
                            ylabel="fraction",xlim=(0,2*pi),size=(400,400*2/3),textsize=4)
                push!(H, h)
                Plot.save((;pyrint=s.pyrint[1], area=s.area[1], cell=s.unit[1])) 
            end
            plot(H[1:min(49,length(H))]..., size=(2000,2000), right_margin = 12Plots.mm)

            Plot.save((;area=sp_area.area[1])) 
        end
    end
end # end dataset loop
