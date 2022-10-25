module manifold

    using Serialization
    using DrWatson
    import Utils
    using Infiltrator
    
    export get_dim_subset
    function get_dim_subset(dict_of_embeds::Dict, dim::Int)::Dict
        Dict(k=>v for (k,v) in dict_of_embeds
             if k[2] == dim)
    end

    function prep_steps()
    end

    export desc_manis
    function desc_manis(;filt, feature_engineer, distance, delim="_")
        if filt !== nothing
            filtstr = "filt=$filt"
        else
            filtstr = "filt=nothing"
        end
        festr   = feature_engineer === nothing ? "feature=nothing" : "feature=$feature_engineer"
        diststr = distance === nothing ? "distance=euclidean" : lowercase("distance=$distance")
        "$(filtstr)$(delim)$(diststr)$(delim)$(festr)"
    end

    export path_manis
    function path_manis(;kws...)
        desc = desc_manis(;kws...)
        datadir("manifold","ca1pfc_manifolds_$desc.serial")
    end

    export save_manis
    function save_manis(;embedding, qlim=[0.02, 0.96], filt, 
                         feature_engineer, distance, data...)
        # Filter
        inds = :inds in keys(data) ? data.inds : Dict()
        for key in keys(embedding)
            inds[key] = Utils.clean.inds_quantile_filter_dims(
                                embedding[key], qlim)
        end
        desc_vars = (;feature_engineer, distance, filt)
        desc      = desc_manis(;desc_vars...)
        data = (;embedding, data..., feature_engineer, distance, filt, inds, desc)
        savefile = path_manis(;feature_engineer, distance, filt)
        @infiltrate
        serialize(savefile, data)
    end

    export load_manis
    function load_manis(mod; feature_engineer, distance, filt)
        savefile = path_manis(;feature_engineer, distance, filt)
        data     = deserialize(savefile)
        for (k, v) in zip(keys(data), values(data))
            Core.eval(mod, :($(Symbol(k)) = $v))
        end
        data
    end


end

