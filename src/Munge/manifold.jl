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
    import Utils
    using Infiltrator
    import Table: CItype
    import Table
    using DataFrames

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
    function load_manis(mod; feature_engineer, distance=:many, filt, tag="")
        savefile = path_manis(;feature_engineer, distance, filt, tag)
        data     = deserialize(savefile)
        skipped=[]
        for (k, v) in zip(keys(data), values(data))
            try
            Core.eval(mod, :($(Symbol(k)) = $v))
            catch
                push!(skipped, k)
            end
        end
        @warn "skipped" skipped
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
            concat::CItype=[])

    Simplifies the process of exploring matched/paired dataframes
    of embeddings annotated/augmented with some other data source
    """
    function EmbeddingFrameFetch(df::DataFrame, pairing::CItype,
            augment::DataFrame=DataFrame(), props::CItype=[];
            concat::CItype=[])
        pairing = pairing isa Vector ? pairing : [pairing]
        if concat != []
            concat = _remainderdims(df, concat)
            @info "concat" concat
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
        @info groupings pairing
        EmbeddingFrameFetch(groupby(df, groupings), pairing, augment, props)
    end
    function Base.getindex(em::EmbeddingFrameFetch, state)
        dfs = em.df[state] 
        [_data(df, em) for df in eachrow(dfs)]
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
    Base.iterate(em::EmbeddingFrameFetch) = em[1], 1
    Base.iterate(em::EmbeddingFrameFetch, state) = count == length(em) ? (em[count+1], count+1) : nothing

end
