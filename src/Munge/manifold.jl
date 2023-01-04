"""
    manifold

module for munging my manifold related analyses

# Major exports

load_manis :: loads manifold data
save_manis :: saves manifold data
path_manis :: gets path to manifold file
desc_manis :: describes manifold data with hyperparmaeter string

"""
module manifold

    using Serialization
    using DrWatson
    using Infiltrator
    using DataFrames, DataFramesMeta

    import Utils
    import Table: CItype
    import Table
    import Filt
    import Utils.namedtup: ntopt_string
    import Utils.dict: load_dict_to_module!

    export get_dim_subset
    function get_dim_subset(dict_of_embeds::Dict, dim::Int)::Dict
        Dict(k=>v for (k,v) in dict_of_embeds
             if k[2] == dim)
    end

    function prep_steps()
    end

    export desc_manis
    function desc_manis(;filt=nothing, 
            feature_engineer=:many, distance=:many, delim="_", tag="")
        if filt !== nothing
            filtstr = "filt=$filt"
        else
            filtstr = "filt=nothing"
        end
        tag = tag == "" ? tag : "$(delim)$tag"
        festr  = feature_engineer === nothing ? "feature=nothing" : "feature=$feature_engineer"
        diststr = distance === nothing ? "distance=euclidean" : lowercase("distance=$distance")
        "$(filtstr)$(delim)$(diststr)$(delim)$(festr)$(tag)"
    end

    export path_manis
    function path_manis(;kws...)
        desc = desc_manis(;kws...)
        datadir("manifold","ca1pfc_manifolds_$desc.serial")
    end

    export save_manis
    function save_manis(;embedding, qlim=[0.02, 0.96], filt, 
                         feature_engineer, distance=:many, tag="", data...)
        # Filter
        inds = :inds in keys(data) ? data.inds : Dict()
        for key in keys(embedding)
            inds[key] = Utils.clean.inds_quantile_filter_dims(
                                embedding[key], qlim)
        end
        desc_vars = (;feature_engineer, distance, filt)
        desc      = desc_manis(;desc_vars...)
        data = (;embedding, data..., feature_engineer, 
                distance, filt, inds, desc)
        savefile = path_manis(;feature_engineer, distance, filt, tag)
        serialize(savefile, data)
    end

    export load_manis
    function load_manis(;feature_engineer, distance=:many, filt, tag="")
        savefile = path_manis(;feature_engineer, distance, filt, tag)
        data     = deserialize(savefile)
        data
    end
    function load_manis(mod; feature_engineer, distance=:many, filt, tag="")
        data = load_manis(;feature_engineer, distance=:many, filt, tag="")
        load_dict_to_module!(mod, data)
        data
    end

    export load_manis_df
    """
        load_manis_workspace

    Handles what Umap_deserialize script used to. Namely, it grabs all of
    the variables used for manifold analyses
    """
    function load_manis_workspace(animal::String, day; 
        filt, areas, distance, feature_engineer, N, trans=:Matrix, kws...)

        data = load_manis(;feature_engineer, filt, distance, tag="$(animal)$(day).$(N)seg")
        embedding_overall = merge(embedding_overall, data[:embedding])
        animal, day = data[:animal], data[:day]
        @time global spikes, beh, ripples, cells = Load.load(animal,day)

        # Filter
        filters = Filt.get_filters()
        if Symbol(filt) in keys(filters)
            beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes; filter_skipmissingcols=true, filters=filters[Symbol(filt)])
        end

        embedding = embedding_overall

        transform(em) = if trans == :DimArray
            (em=em';DimArray(em, (Dim{:time}(beh.time[1:size(em,1)]), Dim{:comp}(1:size(em,2)))))
        elseif trans == :Dataset
            (em=em';Dataset(eachcol(em)...))
        else
            em'
        end

        # Which core would you like to work on?
        #min_dist, n_neighbors, metric, dim = [0.3], [5,150], [:CityBlock, :Euclidean], 3
        min_dist, n_neighbors, metric, dim, feature = [0.3], [150], [:CityBlock], 3, :zscore
        K = filter(k->k.min_dist ∈ min_dist && k.n_neighbors ∈ n_neighbors && k.metric ∈ metric && k.dim == dim && k.feature == feature, keys(embedding))
        embedding = Dict(k=>transform(embedding[k]) for k in K)
        K
        em = make_embedding_df(embedding, inds_of_t, scores, beh)

        Dict(pairs((;K,em)))
    end
    function load_manis_workspace(mod::Module, animal::String, day;kws...)
        data=load_manis_workspace(animal,day;kws...)
        load_dict_to_module!(data)
        data
    end

    export make_embedding_df
    function make_embedding_df(embedding::Dict, inds_of_t::Vector, 
            score::Dict, beh::DataFrame; vars=[])::DataFrame

        em = Table.to_dataframe(embedding, explode=false)
        sc = Table.to_dataframe(score)

        alignkeys = [k for k in (:animal, :day, :n_neighbors, :feature, :metric, :filt, :s, :area, :dim)
                     if k ∈ propertynames(em) && k ∈ propertynames(sc)]
        E, S = groupby(em, alignkeys), 
                 groupby(sc, alignkeys)
        for key in keys(E)
            key = NamedTuple(key)
            e, s = E[key], S[key]
            e.score = s.value
        end

        em[!, :inds_of_t] = inds_of_t[em[!,:s]]
        em[!, :time] = getindex.([beh], em[!,:inds_of_t],[:time])
        for var in vars
            em[!, var] = getindex.([beh], em[!,:inds_of_t], [var])
        end

        em
    end

    export EmbeddingFrameFetch
    mutable struct EmbeddingFrameFetch
        df::GroupedDataFrame
        pairing::CItype
        augment::DataFrame
        props::CItype
    end
    """
        EmbeddingFrameFetch(df::DataFrame, pairing::CItype,
            augment::DataFrame=DataFrame(), props::CItype=[];
            concat::CItype=[],
            ordering::AbstractDict=Dict())

    Simplifies the process of exploring matched/paired dataframes
    of embeddings annotated/augmented with some other data source

    # Inputs
    
    `df` -- embedding dataframe (each embedding is a sample, in the value column,
                    properties of the embedding are metadata about the embedding)
    `pairing` -- column to include multiple frame snippets based on
    `augment` -- dataframe to augment ie pull columns from into the manifold embedding
    `props` -- which properties to pull from augment

    # Optional
    `concat` --- (optional; default=[]) this tells the object to concat over
                 the metadata column `concat` in df
     `ordering` -- (option; default=Dict()) this can instruct the system to first
                    order the grouping keys in the pairing step
    """
    function EmbeddingFrameFetch(df::DataFrame, 
            pairing::CItype, 
            augment::DataFrame=DataFrame(), 
            props::CItype=Not(:); 
            concat::CItype=Not(:), ordering::AbstractDict=Dict())
        pairing = pairing isa Vector ? pairing : [pairing]
        if !isempty(ordering)
            for (sortprop,sortorder) in ordering
                sortprop ∈ propertynames(df) ? 
                sort!(df, sortprop, by=x-> findfirst(x .== sortorder)) : 
                nothing
            end
        end
        if concat != []
            concat = _remainderdims(df, concat)
            @info "EmbeddingFrameFetch" concat
            G = groupby(df, concat)
            K, G′  = keys(G), []
            for k in K
                g = G[k]
                value, score, time, inds_of_t = vcat(g.value...), vcat(g.score...), 
                vcat(g.time...), vcat(collect.(g.inds_of_t)...)
                row = DataFrame(g[1,Not([:value,:score])]) # issue TODO
                row.value, row.score, row.time, row.inds_of_t = [value], [score], 
                                                    [time], [inds_of_t]
                sort!(row,:time)
                push!(G′, row)
            end
            df = vcat(G′...)
        end
        groupings = _remainderdims(df, pairing)
        unique!(groupings)
        @info "EmbeddingFrameFetch" groupings pairing
        EmbeddingFrameFetch(groupby(df, groupings), pairing, augment, 
                            props isa Vector ? union([:time], Symbol.(props)) :
                            props
                           )
    end
    function Base.getindex(em::EmbeddingFrameFetch, state::Int64)
        dfs = em.df[state] 
        [_data(df, em) for df in eachrow(dfs)]
    end
    function Base.getindex(em::EmbeddingFrameFetch, states::UnitRange)
        [em[state] for state in states]
    end
    function _remainderdims(df, pairing)
        pairing = pairing isa Vector ? pairing : [pairing]
        groupings = setdiff(propertynames(df), pairing)
        setdiff(groupings, [:score, :value, :time, :inds_of_t])
    end
    function _data(df::DataFrameRow, em)
        value = df.value
        time  = df.time
        inds_of_t = df.inds_of_t
        data = em.augment[inds_of_t, em.props]
        data.data = collect(eachrow(value))
        data
    end
    Base.length(em::EmbeddingFrameFetch) = length(em.df)
    Base.lastindex(em::EmbeddingFrameFetch) = length(em.df)
    Base.firstindex(em::EmbeddingFrameFetch) = 1
    Base.iterate(em::EmbeddingFrameFetch) = em[1], 1
    Base.iterate(em::EmbeddingFrameFetch, count) = count != length(em) ? 
                                            (em[count+1], count+1) : nothing



end
