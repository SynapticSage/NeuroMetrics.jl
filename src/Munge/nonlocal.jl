module nonlocal

    using DimensionalData, ProgressMeter, DataFrames, DataFramesMeta,
          Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
    using Plots, LazySets, Infiltrator
    using DataStructures: OrderedDict
    
    using DIutils.namedtup
    import DIutils, Random
    import DI: Labels
    import DIutils: Table
    using ...Field: metrics
    using ...Timeshift
    using ...Timeshift.types
    using ...Timeshift.shiftmetrics

    unfiltered_behavior = nothing
    current_metadata = nothing
    function setunfilteredbeh(beh; mod=nonlocal, metadata...) 
        @eval mod unfiltered_behavior = $beh
        @eval mod current_metadata = $metadata
    end
    cuemem_labels = merge(Labels.cuemem, Dict(missing=>"missing"))
    setclab(clabval) = @eval nonlocal cuemem_labels = $clabval
    cueha_labels = merge([Dict(
        (key, 'H') => "H" * replace(value, "cue"=>"C", "mem"=>"M", "nontask"=>"-"),
    (key, 'A') => "A" * replace(value, "cue"=>"C", "mem"=>"M", "nontask"=>"-"))
        for (key, value) in cuemem_labels]...)

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
                    for i in 1:length(U[:hullseg_grid])][DIutils.na, :]
            #hullT = [VPolygon(trans(U[:hullseg_grid][i]))
            #        for i in 1:length(U[:hullseg_grid])][DIutils.na, :]
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
    function get_isolation_summary(spikes, split=[:cuemem]; 
            unfiltered_behavior=unfiltered_behavior, minvel=2.5,
            dropmiss=true)
        if unfiltered_behavior === nothing 
            @error("\nMust either pass in unfiltered behavior or\n"*
                 "set it in this module using `setunfilteredbeh`")
        end
        split = unique([:area, split...])
        iso_sum = combine(
        groupby(dropmissing(spikes,[:isolated,:nearestcyc, :meancyc]), split), 
                          [:isolated,:nearestcyc,:meancyc,:velVec] .=> mean, 
                :velVec=>(x->mean(abs.(x)))=>:speed,
                 (x->nrow(x)))
        if :period ∈ split
            # Calculate time animal spends in each cuemem segment
            task_pers = Table.get_periods(nonlocal.unfiltered_behavior, [:period, :cuemem], 
                                          timefract=:velVec => x->abs(x) > minvel)
            # Total that time and register that column to the isolation summary
            task_pers = combine(groupby(task_pers, [:period, :cuemem]),
                                [:δ,:frac] => ((x,y)->sum(x.*y)) => :timespent,
                                :frac => length => :nsamp)
            DIutils.filtreg.register(task_pers, iso_sum, on="period", transfer=["timespent"])
        elseif :traj ∈ split
            # Calculate time animal spends in each cuemem segment
            task_pers = Table.get_periods(nonlocal.unfiltered_behavior, [:traj, :cuemem], 
                                          timefract=:velVec => x->abs(x) > minvel)
            # Total that time and register that column to the isolation summary
            task_pers = combine(groupby(task_pers, [:period, :cuemem]),
                                [:δ,:frac] => ((x,y)->sum(x.*y)) => :timespent,
                                :frac => length => :nsamp)
            DIutils.filtreg.register(task_pers, iso_sum, on="traj", transfer=["timespent"])
        else
            @info "traj and period not in split"
                                # Calculate time animal spends in each cuemem segment
            task_pers = Table.get_periods(nonlocal.unfiltered_behavior, [:period, :cuemem], 
                              timefract=:velVec => x->abs(x) > minvel)
            dropmissing!(task_pers, :cuemem)
            # Total that time and register that column to the isolation summary
            task_pers = combine(groupby(task_pers, [:cuemem]), 
                                [:δ,:frac] => ((x,y)->sum(x.*y)) => :timespent,
                                :frac => length => :nsamp)
            DIutils.filtreg.register(task_pers, iso_sum, on="cuemem", transfer=["timespent"])
        end

        isempty(iso_sum) ? @error("iso_sum is innapropriately empty") : nothing

        if dropmiss && :cuemem in propertynames(iso_sum)
            dropmissing!(iso_sum, :cuemem)
        end
        if dropmiss && :ha in propertynames(iso_sum)
            dropmissing!(iso_sum, :ha)
        end

        # Add a string column that summarizes cuemem and correct
        if :correct ∈ propertynames(iso_sum) && :cuemem ∈ propertynames(iso_sum)
            iso_sum.cortsk = Vector{Union{String,Missing}}(missing, nrow(iso_sum))
            for (i, (cuemem, correct)) in enumerate(zip(iso_sum.cuemem, iso_sum.correct))
                key = [Int64(cuemem), Int64(correct)]
                if !haskey(Labels.cortsk, key)
                    continue
                end
                inds = findall(x->x.cuemem  == cuemem && 
                                  x.correct == correct, eachrow(iso_sum))
                iso_sum.cortsk[inds] .= Labels.cortsk[key]
            end
        end

        # Acqruire events per time as events  / time spent
        iso_sum = DataFrames.transform(iso_sum, :x1 => :events)[:,Not(:x1)]
        iso_sum.events_per_time = iso_sum.events ./ (iso_sum.timespent)
        iso_sum.cuearea = iso_sum.area .* "\n" .* getindex.([cuemem_labels], iso_sum.cuemem)
        iso_sum.cmlab = getindex.([cuemem_labels], iso_sum.cuemem)
        iso_sum.isolated_events_per_time = iso_sum.isolated_mean .* iso_sum.events_per_time
        iso_sum.iso_event_ratio = iso_sum.isolated_events_per_time ./ iso_sum.events_per_time
        iso_sum.event_iso_ratio = iso_sum.events_per_time ./ iso_sum.isolated_events_per_time 

        if :ha in propertynames(iso_sum)
            K = [iso_sum.cuemem iso_sum.ha] |> eachrow .|> Tuple
            iso_sum.areacuememha = replace.(iso_sum.area .* "\n" .* getindex.([cueha_labels], K),
                                            " "=>"\n")
            iso_sum.cuememha = getindex.([cueha_labels], K)
            iso_sum.cuememhaind = findfirst.(map(x->x.==unique(iso_sum.cuememha), iso_sum.cuememha))
        end


        sortprops = :ha in propertynames(iso_sum) ? [:ha] : []
        ord = Dict("nontask"=>1,"cue"=>2,"mem"=>3, "missing"=>4)
        sort(iso_sum, [sortprops..., :cuearea,DataFrames.order(:cmlab, by=x->ord[x])])
    end

end
