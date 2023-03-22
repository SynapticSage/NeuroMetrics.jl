include("./imports_isolated.jl")
# include(scriptsdir("isolated","imports_isolated.jl"))

datasets = (

            ("RY16",36,:ca1ref), ("RY22",21,:ca1ref),  #("super", 0, :ca1ref),
            # ("RY16",36,:default),("RY22",21,:default), #("super", 0, :default)
           )
# (animal, day, tet) = datasets[2]
@info datasets

opt = isdefined(Main,:opt) ? Main.opt : Dict()
opt = merge(opt, Dict(
           :process_outoffield => false, # process out of field place fields?
           :ploton=>true,
           :matchprops=>[:x,:y,:speed,:startWell,:stopWell],
))

# Plotting turned on?
if opt[:ploton] 
    Plot.on()
else
    Plot.off()
end

init = 1
if init != 1; @warn("initial dataset is $init"); end
(animal, day, tet) = datasets[1]
@info "init dataset" init animal day tet

@showprogress "datasets" for (animal, day, tet) in datasets[init:end]
    @info "loop" animal day tet
    opt["animal"], opt["day"], opt["tet"] = animal, day, tet
    # animal, day, tet = opt["animal"], opt["day"], opt["tet"]
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
    @info "Loading $animal $day"
    @time spikes, beh, cells = DI.load(animal, day, data_source=["spikes","behavior", "cells"])
    GC.gc()
    beh, spikes = DIutils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], 
    filters=filt[datacut], filter_skipmissingcols=true)
    allspikes = copy(spikes)
    beh2 = DI.load_behavior(animal,day)
    Munge.nonlocal.setunfilteredbeh(beh2)
    lfp = DI.load_lfp(animal, day, tet=tet);

    # =============
    ## MUNGE LFP ##
    # =============
    # Center and annotate
    @info "Munging $animal $day"
    begin
        lfp.time = lfp.time .- DI.min_time_records[end]
        # TODO potential bug, 1st time runs, cuts trough-to-trough, second peak-to-peak
        lfp = Munge.lfp.annotate_cycles(lfp, method="peak-to-peak") 
    end
    # Visualize our annotations
    begin
        @assert length(unique(lfp.cycle)) > 1
        #sp = @subset(spikes, :tetrode .== 6);
        cycleplot(lfp)
    end

    # Transfer lfp phase to spikes (for phase locking measurements)
    begin
        DIutils.filtreg.register(lfp, spikes, on="time", transfer=["phase"])
        @assert !all(ismissing.(spikes.phase))
    end
    ##


    # ===================
    # ISOLATED SPIKING
    # ===================
    #Munge.spiking.isolated(sp, lfp, include_samples=true)
    @info "Isospikes $animal $day"
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


    @info "Cycles $animal $day"
    @assert length(unique(lfp.cycle)) > 1
    cycles = Munge.lfp.get_cycle_table(lfp)
    Munge.lfp.annotate_cycles!(spikes, cycles)

    # ENSURING NO SLEEP :: Usually taken care of
    begin
        tsk = DI.load_task(animal, day)
        tsk = DataFrame(unique(eachrow(tsk[:,[:start,:end,:task]])))
        tsk = subset(tsk, :task=>t->t .!= "sleep")
        ismissing.(spikes.velVec)
    end

    # SAVE THE DATA
    begin
        @info "Saving $animal $day"
        filename = path_iso(animal, day, tet)
        !isdir(dirname(filename)) ? mkpath(dirname(filename)) : nothing
        @save "$filename"  animal day lfp spikes allspikes tsk cells beh cycles
    end

    Plot.setparentfolder("isolated")

    # Testing isolation
    Plot.setfolder("phase_locking")
    # All spikes, no cell grouping
    begin
        spikes = :phase ∈ propertynames(spikes) ? spikes[!,Not(:phase)] : spikes
        DIutils.filtreg.register(lfp, spikes, on="time", transfer=["phase"])
        SP = subset(dropmissing(spikes,:isolated),:velVec => v->abs.(v) .> 3)
        SP.pyrint = replace(SP.meanrate .> 4, 0=>"pyr",1=>"int")
        P=[]
        for sp in groupby(SP, :area)
            s = subset(sp, :pyrint => s->s.=="pyr")
            p1=histogram(s.phase, group=s.isolated, normalize=:pdf, alpha=0.5, 
                      legend_title=:isolated, ylabel="fraction",
                      title="PL of $(sp.area[1]) to $tet tetrode"
                     )
            lfp.phase_digitize = Int8.(DIutils.binning.digitize(lfp.phase, 80+1))
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

    begin
        Plot.setfolder("phase_locking", "filtsmooth_$animal-$day-$tet")
        cycles.time = cycles.start
        @info "Registering cycles"
        DIutils.filtreg.register(cycles, spikes; on="time", transfer=["cycle"], match=:prev)
        DIutils.filtreg.register(cycles, lfp; on="time",    transfer=["cycle"], match=:prev)
        @info "Bandstopping broadraw"
        lfp = Munge.lfp.bandstop(lfp, [57, 120], [63, 750]; field=:broadraw, scale=:raw, rounds=1, order=10)
        @info "Smoothing broadraw"
        lfp=Munge.lfp.smooth(lfp, :broadraw; ker=7.5)
        @info "Plotting cycles"
        cycleplot(lfp; otherfield=[:broadraw, :filtbroadraw],
                       kws=[(;linewidth=0.25, linestyle=:dash), 
                            (;linewidth=0.5, linestyle=:dash, legend=:outerbottomright, background_color=:gray)])
        cycleplot(lfp; otherfield=[:broadraw, :smoothbroadraw, :filtbroadraw],
                       kws=[(;linewidth=0.25, linestyle=:dash), (;), (;linewidth=0.5, linestyle=:dash, legend=:outerbottomright, background_color=:gray)])
        Plot.save("filtsmooth")
    end

    # Obtain the firing rate matrix
    R = Munge.spiking.torate(allspikes, beh, gaussian=0.033)

    # Get PFC_units
    pfc_units = @subset(cells,:area.=="PFC").unit
    pfc_units = intersect(pfc_units, collect(R.dims[2]))


    # DATAFRAME-IZE firing rate matrix
    # and add properties of interest
    # --------------------------------
    # Get dataframe of R
    val = :value
    Rdf = DataFrame(R; name=val)
    # Set cycle via start cycle times
    cycles.time = cycles.start;
    cycles.cycle = 1:size(cycles,1);
    DIutils.filtreg.register(cycles, Rdf, on="time", 
                             transfer=["cycle"], 
                             match=:prev)
    matchprops = opt[:matchprops]
    DIutils.filtreg.register(beh, Rdf, on="time", transfer=String.(matchprops))

    commit_vars()


end # end dataset loop
