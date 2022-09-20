module nonlocal

    using Timeshift
    using Timeshift.types
    using Timeshift.shiftmetrics
    using Field.metrics
    using Utils.namedtup
    #using ..Munge.timeshift: getshift
    using DataStructures: OrderedDict
    using DimensionalData
    using ProgressMeter
    using DataFrames, DataFramesMeta
    using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
    using Plots
    using LazySets

    function annotate_nonlocal_spikes!(spikes::DataFrame, unitshift::DimArray, shift::Union{<:Real, Symbol}=0;
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
            cellspace = OrderedDict(r.time => Float32.(element(Singleton(r.x,r.y)))
                            for r in eachrow(spu))
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
            V = collect(values(cellspace))
            infield = V .∈ hull[:,1:hullmax]
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

end
