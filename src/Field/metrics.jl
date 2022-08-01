module metrics

    using StatsBase
    using Entropies: Probabilities
    using NaNStatistics
    using DataFrames
    import DataStructures: OrderedDict
    import ..Field
    import ..Field: ReceptiveField, Grid
    import Table
    import Utils
    import Load.utils: register
    import Load: keep_overlapping_times
    using Images, ImageSegmentation, LazySets
    using Infiltrator

    export Metrics
    export bitsperspike, bitspersecond, coherence, totalcount, maxrate,
           maxcount 

    mutable struct Metrics
        data::AbstractDict{Symbol, Any}
        Metrics()  = new(Dict{Symbol,Any}())
        Metrics(x) = new(x)
    end
    function Base.getindex(M::Metrics, index...)
        Base.getindex(M.data, index...)
    end
    function Base.setindex!(M::Metrics, val, index::Symbol)
        M.data[index] = val
    end
    Base.push!(M::Metrics, p::Pair{Symbol, <:Any}) = push!(M.data, p)
    Base.pop!(M::Metrics, key)::Any = pop!(M.data, key)
    Base.iterate(M::Metrics) = iterate(M.data)
    Base.iterate(M::Metrics, i::Int64) = iterate(M.data, i)
    Base.length(M::Metrics)  = length(M.data)
    function Base.string(S::T where T<:Metrics; sigdigits=2)
        M = ["$k=$(round(v;sigdigits))" for (k,v) in S]
        join(M, " ")
    end

    metric_ban = [:hullzone, :hullsegsizes, :hullseg_inds, :hullseg_grid,
                  :hullseg_inds_cent, :hullseg_grid_cent, :convexhull]

    function Table.to_dataframe(M::Metrics; kws...) 
        kws = (;kws..., key_name=["metric"])
        D = Dict(k=>v for (k,v) in M.data if k ∉ metric_ban || typeof(v) <: AbstractDict)
        Table.to_dataframe(D; kws...)
    end

    function push_metric!(R::ReceptiveField, F::Function; 
            name::Union{Symbol,Nothing}=nothing)
        name = name === nothing ? Symbol(F) : name
        push!(R.metrics, name => F(R))
    end
    pop_metric!(R::ReceptiveField, name::Symbol) = pop!(R.metrics, name)
    function push_metric!(R::ReceptiveField, metrics::AbstractDict)
        for (k, v) in metrics
            push!(R.metrics, k => v)
        end
    end

    skipnan(x) = Iterators.filter(!isnan, x)
    function _convert_to_prob(behProb::AbstractArray)
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
    end

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
    function coherence(F::ReceptiveField)
        rate = F.rate
        rate = (rate .- nanmean(rate))./nanstd(rate)
        neighbors = Vector{Float32}(undef, length(vec(rate)))
        subjects     = Vector{Float32}(undef, length(vec(rate)))
        interaction_mat = get_interaction_mat(ndims(rate))
        interaction_mat = get_interaction_mat(ndims(rate))
        for (i,ind) in enumerate(CartesianIndices(rate))
            subject = rate[ind]
            interactions = collect(Tuple(ind))[Utils.na,:] .+ interaction_mat
            interactions[:,1] = max.(1, min.(size(rate,1), interactions[:,1]))
            interactions[:,2] = max.(1, min.(size(rate,2), interactions[:,2]))
            int_samples = [rate[row...] for row in eachrow(interactions)]
            neighbors[i] = nanmean(int_samples)
            subjects[i] = subject
        end
        pearson = min(nancor(subjects, neighbors), 1f0)
        #@info pearson # sometimes, we get values BARELY above 1
        # project pearson to gaussian for sig
        0.5 * log( (1+pearson) / (1-pearson) )
    end

    function get_interaction_mat(n::Int)
        poss = [1,-1]
        all_possible_neighbors = Iterators.product((poss for i in 1:n)...)
        all_possible_neighbors = collect(all_possible_neighbors)
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

        halfmast  = nanquantile(vec(field.rate), thresh)
        bw        = field.rate .> halfmast
        dist      = 1 .- distance_transform(feature_transform(bw))
        markers   = label_components( (!).(dist .< 0))
        segments  = watershed(dist, markers)
        hullzones = bw .* labels_map(segments)

        # Instruction for how to resort by number of pixels
        instruction = OrderedDict(sort([k=>sum(hullzones .== k) for k in 									                1:maximum(hullzones)], by=x->x[2], rev=true))
        
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
            mets[:hullseg_grid][:toptwohull] = to_grid(mets[:hullseg_inds][:toptwohull], 												field.grid)
            
            mets[:hullseg_grid_cent][:toptwohull] = centroid(iszone, field.grid)
            mets[:hullseg_inds_cent][:toptwohull] = centroid(iszone)
        end

        mets
    end


                                                                    
    # ,---.               o                            |              
    # |---',---.,---..   ..,---.,---.,---.    ,---..  ,|--- ,---.,---.
    # |  \ |---'|   ||   |||    |---'`---.    |---' >< |    |    ,---|
    # `   ``---'`---|`---'``    `---'`---'    `---''  ``---'`    `---^
    #               |                                                 
    #                 
    #               o     ,---.     
    #               .,---.|__. ,---.
    #               ||   ||    |   |
    #               ``   '`    `---'
    #               beyond just the fields

    trajdiversity(spikes::AbstractDataFrame, cell::Int) = length(unique(spikes.traj))
    function trajdiversity(spikes::DataFrame, beh::DataFrame)
        _, spikes = register(beh,spikes;transfer=["traj"])
        combine(groupby(spikes, :traj), trajdiversity)
    end


end
