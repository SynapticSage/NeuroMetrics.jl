module metrics

    using StatsBase
    using Entropies: Probabilities
    using NaNStatistics
    using DataFrames, DataFramesMeta
    import DataStructures: OrderedDict
    using DimensionalData
    import ..Field
    import ..Field: ReceptiveField, Grid
    import Table
    import Table: vec_arrayofarrays!
    import Utils
    import Load.utils: register
    import Load: keep_overlapping_times
    using Images, ImageSegmentation, LazySets
    import TextWrap
    using Infiltrator
    using LazyGrids
    using MATLAB
    using ProgressMeter

    export MetricSet
    export bitsperspike, bitspersecond, coherence, totalcount, maxrate,
           maxcount 
    export push_metric!, pop_metric!
    export run_metrics!
    export push_celltable!, push_spiketable!
    export push_dims!

    mat_initialized = false

    metric_ban = [:hullzone, :hullsegsizes, :hullseg_inds, :hullseg_grid,
                  :hullseg_inds_cent, :hullseg_grid_cent, :convexhull]
    function apply_metric_ban(X::T)::T where T<:AbstractDataFrame
        inds = X.metric .∉ [metric_ban]
        X[inds, :]
    end
    function apply_metric_ban(X::T)::T where T <: AbstractDict
        typeof(X)(k=>v for (k,v) in X
                  if k ∉ metric_ban)
    end


    mutable struct MetricSet
        data::AbstractDict{Symbol, Any}
        MetricSet()  = new(Dict{Symbol,Any}())
        MetricSet(x) = new(x)
    end
    Metrics = MetricSet # for now I require this alias, but eventually, I need to transition my field objects to MetricSet name -- Metrics prevents me from importing Metrics package (third-party)
    Base.getindex(M::MetricSet, index...)            = Base.getindex(M.data, index...)
    Base.setindex!(M::MetricSet, val, index::Symbol) = M.data[index] = val
    Base.push!(M::MetricSet, p::Pair{Symbol, <:Any}) = push!(M.data, p)
    Base.pop!(M::MetricSet, key)::Any                = pop!(M.data, key)
    Base.iterate(M::MetricSet)                       = iterate(M.data)
    Base.iterate(M::MetricSet, i::Int64)             = iterate(M.data, i)
    Base.length(M::MetricSet)                        = length(M.data)
    Base.keys(M::MetricSet)                          = keys(M.data)
    Base.values(M::MetricSet)                        = values(M.data)
    Base.filter(f::Function, M::MetricSet)           = filter(f, M.data)
    Base.getindex(R::ReceptiveField, ind::Symbol)       = Base.getindex(R.metrics, ind)
    Base.setindex!(R::ReceptiveField, val, ind::Symbol) = Base.setindex!(R.metrics, val, ind)
    #Base.setindex(R::ReceptiveField, val, ind)  = Base.setindex(R.metrics, val, ind)
    Base.keys(R::ReceptiveField)     = keys(R.metrics.data)
    Base.values(R::ReceptiveField)   = values(R.metrics.data)
    Base.getindex(R::T where T <: AbstractDimArray{<:ReceptiveField}, 
                  ind::Symbol) = Base.getindex.(R, ind)
    Base.getindex(R::T where T <: AbstractArray{<:ReceptiveField}, 
                  ind::Symbol) = Base.getindex.(R, ind)
    function Base.setindex!(R::T where T <: AbstractArray{<:ReceptiveField}, 
                  val::Union{Real,AbstractArray}, ind::Symbol)
        setindex!.(R, val, ind)
    end
    function Base.setindex!(R::T where T <: AbstractDimArray{<:ReceptiveField}, 
                  val::Union{Real,AbstractArray}, ind::Symbol) 
        R = [r.data for r in R]
        setindex!.(R, val, ind)
    end

    function Base.string(S::T where T<:MetricSet; sigdigits=2, width=40)
        M1 = ["$k=$(round.(v;sigdigits))" for (k,v) in S
             if S isa AbstractArray || typeof(v) <: Real]
        M2 = ["$k=$v" for (k,v) in S
             if  typeof(v) <: AbstractString || typeof(v) <: Symbol]
        TextWrap.wrap(join([M2;M1], " "); width)
    end

    function Table.to_dataframe(M::MetricSet; kws...) 
        kws = (;kws..., key_name=["metric"])
        D = Dict(k=>v for (k,v) in M.data 
                 if k ∉ metric_ban || typeof(v) <: AbstractDict)
        Table.to_dataframe(D; explode=false, kws...)
    end

    """
    unstackMetricDF

    splays out the metric property into columns
    """
    function unstackMetricDF(df::DataFrame)
        if all([:metric,:shift] .∈ [propertynames(df)])
            df = unstack(df, :shift, :metric, :value, allowduplicates=true)
            sort!(df, :shift)
            vec_arrayofarrays!(df)
        end
        df
     end

     """
     remove a configured set of metrics in the module variable metric_ban 
     """
    function apply_metric_ban(X::MetricSet)::MetricSet
        X.data = typeof(X.data)(k=>v for (k,v) in X
                  if k ∉ metric_ban)
        X
    end


    # ----------------- PUSH AND POP FOR METRICS -----------------
    """
        push_metric!

    apply a function with args and kws to generate a metric with a name
    defauling to Symbol(F)
    """
    function push_metric!(R::ReceptiveField, F::Function, args...; 
            name::Union{Symbol}=Symbol(F), kws...)
        #name = name === nothing ? Symbol(F) : name
        push!(R.metrics, name => F(R, args...; kws...))
    end
    """
        push_metric!

    push an existing dict of metrics into metrics
    """
    function push_metric!(R::ReceptiveField, metrics::AbstractDict)
        for (k, v) in metrics
            push!(R.metrics, k => v)
        end
    end
    """
        push_metric!

    push a list of keys.=>values pairs
    """
    function push_metric!(R::ReceptiveField, 
            keys::Union{Base.KeySet,Vector{Symbol}},
            values::Union{Vector, Base.ValueIterator})
        for (k, v) in zip(keys, values)
            push!(R.metrics, k => v)
        end
    end
    """
        push_metric!

    push a value to a symbol
    """
    push_metric!(R::ReceptiveField, key::Symbol, value) = push!(R.metrics, 
                                                                key => value)
    """
        push_metric!

    push to a set of receptive fields
    """
    function push_metric!(R::AbstractArray{ReceptiveField}, F::Function, args...; 
            name::Union{Symbol,Nothing}=nothing, prog::Bool=false, kws...)
        if prog
            prog_name = name===nothing ? String(Symbol(F)) : string(name)
            @info "Starting $prog_name"
            @showprogress prog_name [push_metric!(r, F, args...; name, kws...) for r in R]
        else
            [push_metric!(r, F, args...; name, kws...) for r in R]
        end
    end
    function push_metric!(R::Array{ReceptiveField}, F::Function, args...; 
            name::Union{Symbol,Nothing}=nothing, kws...)
        if prog
            prog_name = name===nothing ? String(Symbol(F)) : string(name)
            @info "Starting $prog_name"
            @showprogress prog_name [push_metric!(r, F, args...; name, kws...) for r in R]
        else
            [push_metric!(r, F, args...; name, kws...) for r in R]
        end
    end
    function push_metric!(R::Vector{ReceptiveField}, 
            keys::Union{Vector{Symbol},Base.KeySet},
            values::Union{Vector, Base.ValueIterator})
        for (r, k, v) in zip(R, keys, values)
            push_metric!(r, k, v)
        end
    end
    pop_metric!(R::ReceptiveField, name::Symbol) = pop!(R.metrics, name)
    Base.push!(R::ReceptiveField, pos...; kws...) = Base.push!(R, pos...;kws...)

    """
        run_metrics!

    function for mass running metrics on a receptive field or set of fields
    """
    function run_metrics!(R::Union{ReceptiveField,
                                  AbstractArray{ReceptiveField}}, instructions)
        for (func, (pos, kws)) in instructions
            push_metric!(R, func, pos...; kws...)
        end
    end
    # ------------------------------------------------------------

    # ---------------- SPECIAL METHODS FOR ADDING METS -----------
    function push_celltable!(r::ReceptiveField, cells::AbstractDataFrame, 
            columns...; cell=r.metrics[:unit])
        for field in columns
            push_metric!(r, field, @subset(cells, :unit .== cell)[1,field])
        end
    end
    function push_celltable!(R::Array{ReceptiveField}, cells::AbstractDataFrame,
            columns...)
        for r in R
            if :unit ∈ R.metrics
                cell = R.metrics[:unit]
            end
            push_celltable!(r, cells, columns...; cell)
        end
    end

    function push_celltable!(R::DimArray{ReceptiveField},
                             cells::AbstractDataFrame, columns...)
        units = R.dims[findfirst(name.(R.dims) .== :unit)]
        units = units .* ones(Broadcast.broadcast_shape(size(units), size(R)))

        for (r,cell) in zip(R, units)
            push_celltable!(r, cells, columns...; cell)
        end
    end

    function push_dims!(R::DimArray{ReceptiveField})
        dims = ndgrid(collect.(R.dims)...)
        N = name(R.dims)
        for (r, D...) in zip(R, dims...)
            for (n, d) in zip(N,D)
                push_metric!(r, n, d)
            end
        end
    end

    function push_dims!(R::DimArray{ReceptiveField}, values::AbstractVector; dim, metric::Symbol)
        N = dim isa Symbol ? name(R.dims) : (1:ndims(R))
        n = findfirst(N.==dim)
        values = Utils.permutevec(values, n)
        values = values .* ones(size(R))
        for (r, v) in zip(R, values)
            push_metric!(r, metric, v)
        end
    end

    function push_spiketable!(R::ReceptiveField, spikes::AbstractDataFrame, field;
            cell)
        push_metric!(R, field, @subset(cells, :unit .== cell)[1,field])
    end
    # ---------------- SPECIAL METHODS FOR ADDING METS -----------

    skipnan(x) = Iterators.filter(!isnan, x)
    function _convert_to_prob(behProb::AbstractArray)
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
    end


    # --------------- CELL WISE METRICS --------------------------

    information(F::T where T <: ReceptiveField; method=:spatialinformation) = information(F.rate, F.occ.prob; method)
    function information(F::Dict, behProb::Probabilities; method=:spatialinformation)
        I = Dict()
        method = eval(method)
        for i in keys(F)
            I[i] = method(F[i], behProb)
        end
        return I
    end
    information(F::Array, behProb::Probabilities; method=:spatialinformation) = eval(method)(F, behProb)


    function bitsperspike(F::T where T<:ReceptiveField)
        bitsperspike(F.rate, F.occ.prob)
    end
    function bitsperspike(firingrate::AbstractArray, behProb::Probabilities)
        firingrate = vec(firingrate)
        FRoverMeanFR = firingrate ./ nanmean(firingrate)
        # paper r = ∑ pᵢ * rᵢ
        # what I've written here: r = ∑ R_occᵢ / N = μ(R_occᵢ), not the actual rate, but occ norm rate
        I = nansum( behProb .* FRoverMeanFR .* log2.(FRoverMeanFR) )
        # what i've written here: ∑ p(state) * r(state)/r̅ * log₂( r(state) / r̅ )
        # paper : ∑ pᵢ * ( rᵢ/r̅ ) * log₂( rᵢ / r̅ )
        return I
    end
    spatialinformation = bitsperspike
    function bitspersecond(firingrate::AbstractArray, behProb::Probabilities)
        firingrate = vec(firingrate)
        FRoverMeanFR = firingrate ./ nanmean(firingrate)
        # paper r = ∑ pᵢ * rᵢ
        # what I've written here: r = ∑ R_occᵢ / N = μ(R_occᵢ), not the actual rate, but occ norm rate
        I = nansum( behProb .* firingrate .* log2.(FRoverMeanFR) )
        # what i've written here: ∑ p(state) * r(state)/r̅ * log₂( r(state) / r̅ )
        # paper : ∑ pᵢ * ( rᵢ/r̅ ) * log₂( rᵢ / r̅ )
        return I
    end

    function kraskov_mutualinformation(spikeRate, behRate; kws...)
        K = Entropies.Kraskov(;kws...) 
        D = DataSet(spikeRate, behCount)
        HXY = Entropies.genentropy(D, K)
        HX  = Entropies.genentropy(spikeRate, K)
        HY  = Entropies.genentropy(behRate, K)
        HX + HY - HXY
    end
    
    """
    mutualinformation

    get the mutual information of a single cell
    """
    function mutualinformation(beh, spikes, props, unit)
        beh, spikes = keep_overlapping_times(beh, spikes)
    end

    function totalcount(F::ReceptiveField)
        sum(F.count)
    end
    function maxcount(F::ReceptiveField)
        maximum(F.count)
    end
    function maxrate(F::ReceptiveField)
        nanmaximum(F.rate)
    end
    function meanrate(F::ReceptiveField)
        nanmean(F.rate)
    end

    function hull_withimages(X::BitMatrix)::Vector{CartesianIndex}
        convexhull(X)
    end
    function hull_withlazysets(X::BitMatrix)::Vector{Vector{Float32}}
        collect(convex_hull(array_of_arrays(findall(X))))
        #[convert(Vector{Float32}, x) for x in h]
    end
    function centroid(X::BitArray)::Vector{Int32}
        round.(mean(array_of_arrays(findall(X))))
    end
    function centroid(X::BitArray, grid::Grid)::Vector{Float32}
        grid.grid[ centroid(X)... ]
    end
    function centroid(X::ReceptiveField)::Vector{Float32}
        filled = (!).(isnan.(X.rate))
        sum(X.grid.grid[filled] .* X.rate[filled]) ./ sum(X.rate[filled])
    end
    array_of_tuples(X::Vector{<:CartesianIndex}) = [x.I for x in X]
    array_of_arrays(X::Vector{<:CartesianIndex}) = [collect(x.I) for x in X]
    array_of_singleton(X::Vector{<:CartesianIndex}) = [Singleton(collect(x.I)) for x in X]
    
    loopup_coord(c::Tuple, F::Field.ReceptiveField) = F.grid.grid[c...]
    to_grid(X::Vector{<:Union{Tuple,Vector}}, grid::T where T <: Grid) = [grid.grid[Int32.(x)...] for x in X]


    function argmax(X::ReceptiveField)::Vector{Float32}
        filled = (!).(isnan.(X.rate))
        i = Base.argmax(X.rate[filled])
        X.grid.grid[i]
    end


    """
        coherence(F::Field.ReceptiveField)

    computes the spatial coherence of a field
    """
    function jcoherence(F::ReceptiveField;
             zscore::Bool=true, 
             skip_edge::Bool=false,
            convert_to_sig::Bool=true)

        if ndims(F.rate) > 2
            return NaN
        end


        rate = F.rate
        rate = zscore ? (rate .- nanmean(rate))./nanstd(rate) : rate

        neighbors = Vector{Float32}(undef, length(vec(rate)))
        subjects  = Vector{Float32}(undef, length(vec(rate)))

        interaction_mat = get_interaction_mat(ndims(rate))

        for (i,ind) in enumerate(CartesianIndices(rate))
            # Skip the edge of the track
            if skip_edge && is_edge_pixel(ind.I, rate)
                neighbors[i] = subjects[i] = NaN
                continue
            end

            # Get the subject pixel
            subject = rate[ind]

            # OBTAIN INDICES OF INTERACTION
            interactions = collect(Tuple(ind))[Utils.na,:] .+ interaction_mat
            interactions[:,1] = max.(1, min.(size(rate,1), interactions[:,1]))
            interactions[:,2] = max.(1, min.(size(rate,2), interactions[:,2]))

            # PULL NEIGHBORS for a pixel
            int_samples = [rate[row...] for row in eachrow(interactions)]

            # Store mean neighbor and subject
            neighbors[i] = nanmean(int_samples)
            subjects[i]  = subject
        end

        # Throw out any nans (this isn't necessary)
        good_sample = (!).(isnan.(neighbors) .|| isnan.(subjects))
        subjects, neighbors = subjects[good_sample], neighbors[good_sample]

        pearson = min(cor(subjects, neighbors), 1f0)
        convert_to_sig ? fisher(pearson) : pearson
    end

    function coherence(F::ReceptiveField;
             skip_edge::Bool=false,
            convert_to_sig::Bool=true)


        rate, count, occ = F.rate, F.count, F.occ.count
        if ndims(rate) != 2
            return NaN
        end


        neighbors = Vector{Float32}(undef, length(vec(rate)))
        subjects  = Vector{Float32}(undef, length(vec(rate)))

        interaction_mat = get_interaction_mat(ndims(rate))

        for (i,ind) in enumerate(CartesianIndices(rate))
            # Skip the edge of the track
            if skip_edge && is_edge_pixel(ind.I, rate)
                neighbors[i] = subjects[i] = NaN
                continue
            end

            # Get the subject pixel
            subject = rate[ind]

            # OBTAIN INDICES OF INTERACTION
            interactions = collect(Tuple(ind))[Utils.na,:] .+ interaction_mat
            interactions[:,1] = max.(1, min.(size(rate,1), interactions[:,1]))
            interactions[:,2] = max.(1, min.(size(rate,2), interactions[:,2]))

            # PULL NEIGHBORS for a pixel
            neigh_spikes = nansum([count[row...] for row in eachrow(interactions)])
            neigh_occ    = nansum([occ[row...] for row in eachrow(interactions)])

            # Store mean neighbor and subject
            neighbors[i] = neigh_spikes/neigh_occ
            subjects[i]  = subject
        end

        # Throw out any nans (this isn't necessary)
        good_sample = (!).(isnan.(neighbors) .|| isnan.(subjects))
        subjects, neighbors = subjects[good_sample], neighbors[good_sample]

        pearson = min(cor(subjects, neighbors), 1f0)
        convert_to_sig ? fisher(pearson) : pearson
    end


    function fisher(pearson::Real)
       0.5 * log( (1+pearson) / (1-pearson) )
    end

    function is_edge_pixel(ind::Tuple, rate::Array)::Bool
        edge = false
        for coord in ind
            if (coord == 1) || (coord == size(rate,coord))
                edge = true
                break
            end
        end
        edge
    end

    function blake_coherence(F::ReceptiveField; zscore::Bool=true)
        metrics.mat_initialized ? nothing : initialize_mat()
        mat"blake_coherence(double($(F.count)), double($(F.occ.count)), $(size(F.count,1)), $(size(F.count,2)))"
    end
    function jake_coherence_spearman(F::ReceptiveField; significance::Bool=true)
        metrics.mat_initialized ? nothing : initialize_mat()
        val = mat"jake_cohenrence(double($(F.rate)), double($(F.occ.count)), 'Spearman')"
        significance ? fishÿner(val) : val
    end
    function jake_coherence_pearson(F::ReceptiveField; significance::Bool=true)
        metrics.mat_initialized ? nothing : initialize_mat()
        val = mat"jake_coherence(double($(F.rate)),double($(F.occ.count)), 'Pearson')"
        significance ? fisher(val) : val
    end
    function initialize_mat()
        @eval metrics mat_initialized = true
        mat"addpath('~/Code/metrics')"
    end


    function get_interaction_mat(n::Int)
        poss = [1,0,-1]
        all_possible_neighbors = Iterators.product((poss for i in 1:n)...)
        all_possible_neighbors = [x for x in collect(all_possible_neighbors)
                                    if x != (0,0)]
        all_possible_neighbors = vcat(all_possible_neighbors...)
        vcat([collect(x)[Utils.na, :] for x in all_possible_neighbors]...)
    end

    #                                                                  
    #            ,---.                             |         |    |    
    #            |    ,---.,---..    ,,---..  ,    |---..   .|    |    
    #            |    |   ||   | \  / |---' ><     |   ||   ||    |    
    #            `---'`---'`   '  `'  `---''  `    `   '`---'`---'`---'
    #                                                                  

	struct HullSet
		hulls::Dict
	end
	function Base.:∈(H::HullSet, x)
		any([x ∈ v for v in values(H)])
	end
	function Base.:⊆(H::HullSet, x)
		any([x ⊆ v for v in values(H)])
	end
	Base.iterate(H::HullSet) = Base.iterate(H.hulls)
	Base.iterate(H::HullSet, i::Int64) = Base.iterate(H.hulls, i::Int64)
	Base.keys(H::HullSet) = Base.keys(H.hulls)
	function plothullset!(H::HullSet)
		for i in sort(collect(filter(x-> x isa Int, keys(H))))
			plot!(VPolygon(H.hulls[i]))
			loc = Iterators.flatten(mean(H.hulls[i],dims=1))
			annotate!(loc..., text(string(i), :white))
		end
    end

    function resort_zones(hullzones, instruction::OrderedDict)
        newhullzones = copy(hullzones)
        segsizes = Vector{Int32}(undef, length(instruction))
        for (new, (curr, count)) in enumerate(instruction)
            #if current != new
            #	@info "different"
            #end
            hzinds = hullzones .== Int64(curr)
            #@info curr any(hzinds)
            segsizes[new]         = sum(convert(Array{Int32}, hzinds))
            newhullzones[hzinds] .= new
        end
        return newhullzones, 
               OrderedDict(zip(1:length(instruction), 
                           values(instruction))),
               OrderedDict(zip(1:length(instruction),
                               segsizes))
    end

    function convexhull(field::ReceptiveField;
            thresh=0.85, tophull=Inf, toptwohulls::Bool=false)

        if ndims(field.rate) < 2
            return Dict()
        end

        halfmast  = nanquantile(vec(field.rate), thresh)
        bw        = field.rate .> halfmast
        dist      = 1 .- distance_transform(feature_transform(bw))
        markers   = label_components( (!).(dist .< 0))
        segments  = watershed(dist, markers)
        hullzones = bw .* labels_map(segments)

        # Instruction for how to resort by number of pixels
        instruction = OrderedDict(sort([k=>sum(hullzones .== k) for k in 									               
                                        1:maximum(hullzones)], by=x->x[2], rev=true))
        
        newhullzones, resorted, segsizes = resort_zones(hullzones, instruction)
        
        mets = Dict()
        mets[:hullzone]      = newhullzones
        mets[:hullsegsizes]  = segsizes
        
        ordered_seg = sort(collect(keys(segments.segment_pixel_count)))
        zones = 1:maximum(newhullzones)
        mets[:hullsegsizes] = OrderedDict(k=>segments.segment_pixel_count[k]
                                          for k in ordered_seg)
        mets[:hullseg_inds] = Dict{Union{Int, Symbol},Any}()
        mets[:hullseg_grid] = Dict{Union{Int, Symbol},Any}()
        mets[:hullseg_inds_cent] = Dict{Union{Int, Symbol},Any}()
        mets[:hullseg_grid_cent] = Dict{Union{Int, Symbol},Any}()
            
        for zone in filter(zone->zone < tophull, zones)
            iszone = newhullzones .== zone
            mets[:hullseg_inds][zone] = hull_withlazysets(iszone)
            mets[:hullseg_grid][zone] = to_grid(mets[:hullseg_inds][zone], 
                                                field.grid)
            mets[:hullseg_inds_cent][zone] = centroid(iszone)
            mets[:hullseg_grid_cent][zone] = centroid(iszone, field.grid)
        end

        # Combine the top two hulls
        if toptwohulls
            iszone = newhullzones .== ordered_seg[1] .||
                     newhullzones .== ordered_seg[2]

            mets[:hullseg_inds][:toptwohull] = hull_withlazysets(iszone)
            mets[:hullseg_grid][:toptwohull] = to_grid(mets[:hullseg_inds][:toptwohull],
                                                       field.grid)
            
            mets[:hullseg_grid_cent][:toptwohull] = centroid(iszone, field.grid)
            mets[:hullseg_inds_cent][:toptwohull] = centroid(iszone)
        end

        mets
    end
    
    #trajdiversity(spikes::AbstractDataFrame, cell::Int) = length(unique(spikes.traj))
    function trajdiversity(spikes::DataFrame, beh::DataFrame)
        _, spikes = register(beh,spikes;transfer=["traj"])
        combine(groupby(spikes, :traj), trajdiversity)
    end
    trajdiversity(R::ReceptiveField, spikes::AbstractDataFrame;
                  cell=R.metrics[:unit]) = trajdiversity(spikes, cell)


    #function fromcelltable(R::DimArray{ReceptiveField}, cells::AbstractDataFrame,
    #        columns...)
    #    for r in R
    #        if :unit ∈ R.metrics
    #            cell = R.metrics[:unit]
    #        end
    #        push_metric!(r, cells, columns...; cell)
    #    end
    #end

    #mutable struct MetricDF
    #    data::DataFrame
    #    MetricDF()  = new(DataFrame())
    #    MetricDF(x::DataFrame) = new(x)
    #end
    #function Table.to_dataframe(M::MetricDF; kws...) 
    #    kws = (;kws..., key_name=["metric"])
    #    D = M.data[M.data.metric ∉ [metric_ban],:]
    #    Table.to_dataframe(D; explode=false)
    #end

end

#    function push_celltable!(r::ReceptiveField, cells::AbstractDataFrame, 
#            columns...; cell=r.metrics[:unit])
#        for field in columns
#            push_metric!(r, field, @subset(cells, :unit .== cell)[1,field])
#        end
#    end
#    function push_celltable!(R::Array{ReceptiveField}, cells::AbstractDataFrame,
#            columns...)
#        for r in R
#            if :unit ∈ R.metrics
#                cell = R.metrics[:unit]
#            end
#            push_celltable!(r, cells, columns...; cell)
#        end
#    end
#
#    function push_celltable!(R::DimArray{ReceptiveField},
#                             cells::AbstractDataFrame, columns...)
#        @infiltrate
#        units = R.dims[findfirst(name.(R.dims) .== :unit)]
#        units = units .* ones(Broadcast.broadcast_shape(size(units), size(R)))
#
#        for (r,cell) in zip(R, units)
#            push_celltable!(r, cells, columns...; cell)
#        end
#    end
#
#    function push_dims!(R::DimArray{ReceptiveField})
#        dims = ndgrid(collect.(R.dims)...)
#        N = name(R.dims)
#        for (r, D...) in zip(R, dims...)
#            for (n, d) in zip(N,D)
#                push_metric!(r, n, d)
#            end
#        end
#    end
#
#    function push_spiketable!(R::ReceptiveField, spikes::AbstractDataFrame, field;
#            cell)
#        push_metric!(R, field, @subset(cells, :unit .== cell)[1,field])
#    end
#
#
