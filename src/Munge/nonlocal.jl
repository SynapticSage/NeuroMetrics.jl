module nonlocal

    using Timeshift
    using Timeshift.types
    using Timeshift.shiftmetrics
    using Field.metrics
    using Utils.namedtup
    import Utils
    import Table
    using DataStructures: OrderedDict
    using DimensionalData
    using ProgressMeter
    using DataFrames, DataFramesMeta
    using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
    using Plots
    using LazySets
    import Random
    using Infiltrator

    beh2 = nothing
    setunfilteredbeh(beh) = @eval nonlocal beh2 = $beh
    clab = nothing
    setclab(clabval) = @eval nonlocal clab = $clabval

    export annotate_nonlocal_spikes!
    function annotate_nonlocal_spikes!(spikes::DataFrame, cells::DataFrame, unitshift::DimArray, shift::Union{<:Real, Symbol}=0;
                                        thresh=0.8, hullmax=1,
                                        n_examples=20, shuf=true)
        # Push details about fields and best shift
        push_celltable!(unitshift, cells, :unit, :area)
        push_dims!(unitshift)
        push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)
        shift0 = unitshift[shift=At(shift)]
        # Add convex hull
        push_metric!(shift0, metrics.convexhull; thresh=0.80)
        hsg = [x[:convexhull][:hullseg_grid] for x in shift0];
        shift0[:hullseg_grid] = hsg;
        shift0 = DimArray(shift0, Dim{:unit}(shift0[:unit]));
        # Obtain nonlocal spikes
        P = []
        spikes.infield .= Vector{Union{Missing,Bool}}(undef, size(spikes,1))
        #trans(X) = [x[end:-1:begin] for x in X]
        for unit in shift0[:unit]
            inds = findall(spikes.unit .== unit .&& (!).(ismissing.(spikes.x)))
            if isempty(inds) || shift0[unit=At(unit)][:totalcount] == 0; 
            @info "skipping" unit
            continue; 
            else
            @info "processing" unit
            end
            spu = view(spikes, inds, :)
            histogram2d(spu.x, spu.y)
            #epsilon = eps(Float64)*2
            #r.time .+= ones(size(r.time))
            cellspace = [Float32.(element(Singleton(r.x,r.y)))
                         for r in eachrow(spu)] # the transform not guarenteeed same size
            U = shift0[unit=At(unit)]
            hull = [VPolygon(U[:hullseg_grid][i])
                    for i in 1:length(U[:hullseg_grid])][Utils.na, :]
            #hullT = [VPolygon(trans(U[:hullseg_grid][i]))
            #        for i in 1:length(U[:hullseg_grid])][Utils.na, :]
            if isempty(hull)
                @warn "no hull"
                continue
            else
                #@infiltrate
            end
            V = collect(cellspace)
            #@infiltrate unit == 92
            infield = V .∈ hull[:,1:min(size(hull,2),hullmax)]
            #infieldT = V .∈ hullT
            spu.infield = vec(any(infield,dims=2))
            title="unit=$unit,\nmean infield=$(round(mean(infield), sigdigits=2))"
            push!(P, Plots.plot(U; transpose=true))
            for h in 1:hullmax
                Plots.plot!(hull[h]; title)
            end
        end
        if shuf == true 
            Random.shuffle(P) 
        end
        Plots.plot(P[1:n_examples]...; 
        colorbar=false, titlefontsize=6, tickfontsize=6, framestyle=:none, size=(800,800))
    end

    export get_isolation_summary
    """
        get_isolation_summary

    acquires an isolation summary table for a spikes DataFrame

    # Input
    spikes  <- spikes dataframe
    split   <- group/split in groupby to create the summary table

    # Requires
    You already added a column specifying if each spike is isolated

    # Returns
    A summary table
    """
    function get_isolation_summary(spikes,split=[:cuemem])
        iso_sum = combine(groupby(dropmissing(spikes,[:isolated,:nearestcyc, :meancyc]), [:area, split...]), 
                          [:isolated,:nearestcyc,:meancyc,:velVec] .=> mean, (x->nrow(x)))
        @infiltrate
        if :period ∈ split
            # Calculate time animal spends in each cuemem segment
            task_pers = Table.get_periods(beh2, [:period, :cuemem], 
                                          timefract=:velVec => x->abs(x) > 2)
            # Total that time and register that column to the isolation summary
            task_pers = combine(groupby(task_pers, [:period, :cuemem]),
                                [:δ,:frac] => ((x,y)->sum(x.*y)) => :timespent)
            Utils.filtreg.register(task_pers, iso_sum, on="period", transfer=["timespent"])
        elseif :traj ∈ split
            # Calculate time animal spends in each cuemem segment
            task_pers = Table.get_periods(beh2, [:traj, :cuemem], 
                                          timefract=:velVec => x->abs(x) > 2)
            # Total that time and register that column to the isolation summary
            task_pers = combine(groupby(task_pers, [:period, :cuemem]),
                                [:δ,:frac] => ((x,y)->sum(x.*y)) => :timespent)
            Utils.filtreg.register(task_pers, iso_sum, on="traj", transfer=["timespent"])
        else
            # Calculate time animal spends in each cuemem segment
            task_pers = Table.get_periods(beh2, [:period, :cuemem], 
                              timefract=:velVec => x->abs(x) > 2)
            dropmissing!(task_pers, :cuemem)
            # Total that time and register that column to the isolation summary
            task_pers = combine(groupby(task_pers, [:cuemem]), 
                                [:δ,:frac] =>
                                ((x,y)->sum(x.*y)) => :timespent)
            Utils.filtreg.register(task_pers, iso_sum, on="cuemem", transfer=["timespent"])
        end
        # Acqruire events per time as events  / time spent
        iso_sum = transform(iso_sum, :x1 => :events)[:,Not(:x1)]
        iso_sum.events_per_time = iso_sum.events ./ (iso_sum.timespent)
        iso_sum.cuearea = iso_sum.area .* "\n" .* getindex.([clab], iso_sum.cuemem)
        iso_sum.cmlab = getindex.([clab], iso_sum.cuemem)
        iso_sum.isolated_events_per_time = iso_sum.isolated_mean .* iso_sum.events_per_time
        ord = Dict("nontask"=>1,"cue"=>2,"mem"=>3)
        @infiltrate
        sort(iso_sum, [DataFrames.order(:cmlab, by=x->ord[x]),:cuearea])
    end

end
