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
    using ArgParse
    using Statistics
    using MultivariateStats
    import StatsAPI
    using StatsBase
    import DIutils.plotutils: get_lims
    import DIutils.arr: get_quantile_filtered

    import DI
    import DI: Filt
    using DIutils
    using DIutils.Table 
    import DIutils.namedtup: ntopt_string

    """
        parse(args=nothing; return_parser::Bool=false)

    return command line parser for flags that control my manifold
    related analyses
    """
    function parse(args=nothing; return_parser::Bool=false)
        parser = ArgParseSettings()
        default_splits = 10
        default_sps = 10
        default_N = -1
        @add_arg_table parser begin
            "--animal"
                help = "the animal to run default: the super animal"
                arg_type = String
                default = "RY16"
            "--day"
                help = "the day to run"
                arg_type = Int
                default = 36
            "--dataset"
                 help = "dataset preset"
                 arg_type = Int
                 default = 0
            "--splits", "-s"
                help = "splits of the dataset"
                default = default_splits
                arg_type = Int
            "--sps", "-S"
                help = "samples per split"
                default = default_sps
                arg_type = Int
            "--N"
                help = "number of manis"
                default = default_N
                arg_type = Int
            "--filt", "-f"
                help = "type of dataset filter in Filt.jl"
                arg_type = Symbol
                default = :all
            "--feature_engineer", "-F"
                help = "types of features"
                default = :many
                arg_type = Symbol
            "--distance", "-d"
                help = "distance metric(s)"
                default= :many
                arg_type=Symbol
            "--save_frequency",
                help = "how often to save the manis"
                default = 100
                arg_type = Int
        end
        if return_parser
            return parser
        else
            if args !== nothing
                opt = postprocess(parse_args(args, parser))
            else
                opt = postprocess(parse_args(parser))
            end
            @info "parsed options" opt
            return opt
        end
    end

    function postprocess(opt::Dict)
        if opt["N"] == -1
            opt["N"] = opt["sps"] * opt["splits"]
        end
        opt
    end

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
                         feature_engineer, distance=:many, tag, data...)
        # Filter
        inds = :inds in keys(data) ? data.inds : Dict()
        for key in keys(embedding)
            inds[key] = DIutils.clean.inds_quantile_filter_dims(
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
    function load_manis(;feature_engineer, distance=:many, filt, tag)
        @info "load_manis" feature_engineer distance filt tag
        savefile = path_manis(;feature_engineer, distance, filt, tag)
        if isfile(savefile)
            @info "load_manis located savefile" savefile
        else
            @error "load_manis did not locate savefile" savefile tag
        end
        data     = deserialize(savefile)
        data
    end
    function load_manis(mod::Module; feature_engineer, distance=:many, filt, tag)
        data = load_manis(;feature_engineer, distance=:many, filt, tag)
        dict.load_dict_to_module!(mod, Dict(pairs(data)))
        data
    end

    export load_manis_workspace
    """
        load_manis_workspace

    Handles what Umap_deserialize script used to. Namely, it grabs all of
    the variables used for manifold analyses

    Calls load_manis, pulling computing manifolds, indices, time, et cetera. It then
    elaborates this set of variables with a dataframe version of the manifold data,
    individual samples in rows and properties in columns.
    """
    function load_manis_workspace(animal::String, day; 
        filt, distance, feature_engineer, N, trans=:Matrix, embedding_overall=nothing, kws...)

        data = load_manis(;feature_engineer, filt, distance, tag="$(animal)$(day).$(N)seg")

        embedding_overall = embedding_overall === nothing ?
            typeof(data[:embedding])() : embedding_overall
        embedding_overall = merge(embedding_overall, data[:embedding])
        # @time spikes, beh, ripples, cells = DI.load(animal,day)
        beh = DI.load_behavior(animal,day)

        # Filter
        # filters = Filt.get_filters()
        # if Symbol(filt) in keys(filters)
        #     beh, spikes = DIutils.filtreg.filterAndRegister(beh, spikes; filter_skipmissingcols=true, filters=filters[Symbol(filt)])
        # end

        embedding = embedding_overall

        transform(em) = if trans == :DimArray
            (em=em';DimArray(em, 
                             (Dim{:time}(beh.time[1:size(em,1)]),
                              Dim{:comp}(1:size(em,2)))
            ))
        elseif trans == :Dataset
            (em=em';Dataset(eachcol(em)...))
        else
            em'
        end

        # Which core would you like to work on?
        #min_dist, n_neighbors, metric, dim = [0.3], [5,150], [:CityBlock, :Euclidean], 3
        #min_dist, n_neighbors, metric, dim, feature = [0.3], [150], [:CityBlock], 3, :zscore
        #filtkeys = filter(k->k.min_dist ∈ min_dist && k.n_neighbors ∈ n_neighbors && 
        #           k.metric ∈ metric && k.dim == dim && k.feature == feature, 
        #           keys(embedding)) # filter keys
        # embedding = Dict(k=>transform(embedding[k]) for k in filtkeys) # transform embedding entries
        emkeys = keys(embedding)
        embedding = Dict(k=>transform(embedding[k]) for k in emkeys) # transform embedding entries
        emdf = make_embedding_df(embedding, data[:inds_of_t], data[:scores], beh) # get embedding dataframe
        merge(Dict(pairs(data)), Dict(pairs((;emkeys,emdf))))
    end
    function load_manis_workspace(mod::Module, animal::String, day; kws...)
        data=load_manis_workspace(animal,day;kws...)
        dict.load_dict_to_module!(mod, data)
        data
    end

    export make_embedding_df
    """
        make_embedding_df(embedding::Dict, inds_of_t::Vector, score::Dict,
        beh::DataFrame; vars=[], quantile_bounds=(0.001, 0.999))::DataFrame

    # Arguments
    - `embedding::Dict`: Dictionary of embeddings, with keys of the form
        `:animal_day_n_neighbors_feature_metric_filt_s_dataset_dim`
    - `inds_of_t::Vector`: Vector of indices of the time points in the embedding
    - `score::Dict`: Dictionary of scores, with keys of the form
        `:animal_day_n_neighbors_feature_metric_filt_s_dataset_dim`
    - `beh::DataFrame`: DataFrame of behavior, with columns `:animal`, `:day`,
        `:time`, and `:rip_sec`
    - `vars::Vector`: Vector of variables to include in the dataframe. If empty,
    
    # Returns
    - `DataFrame`: Dataframe of the embedding data, with each sample in a row
        and each property in a column. The properties are the keys of the 
        embedding dictionary. The dataframe also has a column for the score 
        of the sample and a column for the behavior of the sample.
    """
    function make_embedding_df(embedding::Dict, inds_of_t::Vector, 
            score::Dict, beh::DataFrame; vars=[],
            quantile_bounds=(0.001, 0.999))::DataFrame

        df = Table.to_dataframe(embedding, explode=false)
        sc = Table.to_dataframe(score)

        alignkeys = [k for k in (:animal, :day, :n_neighbors, :feature, :metric, :filt, :s, :dataset, :dim)
                     if k ∈ propertynames(df) && k ∈ propertynames(sc)]
        df = subset(df, :dataset => a -> (!).(ismissing.(a)))
        assignarea(x::Symbol) = begin
            x = String(x)
            if occursin(x, "ca1")
                :ca1
            elseif occursin(x, "pfc")
                :pfc
            else
                :none
            end
        end
        if :area ∉ propertynames(df) && :dataset ∈ propertynames(df)
            df[!, :area]  = assignarea.(df[!, :dataset])
        end
        E, S = groupby(df, alignkeys), 
               groupby(sc, alignkeys)
        for key in keys(E)
            key = NamedTuple(key)
            try
                e, s = E[key], S[key]
                e.score = s.value
            catch
                @warn "key failed" key
            end
        end

        df[!, :inds_of_t] = inds_of_t[df[!,:s]]
        df[!, :T_partition] = df[!,:s]
        df[!, :T_start], df[!, :T_end] = first.(df[!, :inds_of_t]),
                                         last.(df[!, :inds_of_t])
        df[!, :time] = getindex.([beh], df[!,:inds_of_t],[:time]) # TODO
        for var in vars
            df[!, var] = getindex.([beh], df[!,:inds_of_t], [var])
        end

        # TODO set quantile axis limits
        df.bounds = Vector{Tuple}(undef, size(df,1))
        for row in eachrow(df)
            if !(quantile_bounds isa Vector{Tuple})
                qb = fill(quantile_bounds, row.dim,1)
            else
                qb = quantile_bounds
            end
            row.bounds = ([quantile(v, q) for (v,q) in zip(eachcol(row.value), qb)]...,)
        end

        df
    end

    ########################################################
    ############### Triggered Embedding ####################
    ########################################################

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



    ############################################
    ####### Projection and Embedding ###########
    ############################################
    
    export unity_scale
    """
        unity_scale(B)
    Scale the matrix `B` columns to be between 0 and 1.
    """
    function unity_scale(B)
        B = (B .- minimum(B, dims=1)) ./ (maximum(B, dims=1) .- minimum(B, dims=1))
        # Center 
        B = B .- mean(B, dims=1)
    end

    export get_B_matrix
    """
        get_B_matrix(beh, inds_of_t, s; beh_vars, dummvars)
    Get the behavior matrix for the data fraction `s` from the behavior data `beh`.
    # Arguments
    - `beh`: The behavior data
    - `inds_of_t`: The indices of the data fraction
    - `s`: The data fraction
    - `beh_vars`: The behavior variables to include
    - `dummvars`: The dummy variables to include
    - `weighting` : the weighting of each behavior variable after
            scaling. Default is to weight all variables equally.
            the weighting index across beh_vars and lumped entries
            of each dummary variable.
    # Returns
    - `B`: The behavior matrix
    """
    function get_B_matrix(beh::DataFrame, inds_of_t, s;
    beh_vars = [:x, :y, :cuemem], dummvars=[:stopWell, :startWell], scale=true, 
        scaling=ones(length(beh_vars)+length(dummvars)))
        lb = length(beh_vars)
        B = Matrix(beh[:, beh_vars])
        S = ones(size(B)) .* DIutils.arr.atleast2d(scaling[1:lb])'
        for (v,var) in enumerate(dummvars)
            # Create a dummy code for behavior stopWells
            DS = DIutils.statistic.dummycode(beh[:,var])
            stopWellVars = [Symbol("$(var)_$i") for i in 1:size(DS,2)]
            B = hcat(B, DS)
            sc = ones(size(DS)) .* scaling[lb+v]
            S = hcat(S, sc)
        end
        # Subset for this data fraction
        scale ? 
            unity_scale(B[inds_of_t[s], :]) .* S[inds_of_t[s],:] : 
            B[inds_of_t[s], :]
    end

    export project_onto_behavior
    """
            project_onto_behavior(R::Matrix, B::Matrix)"
    Project the data `R` onto the behavior variables `beh_vars` in `beh`.
    # Arguments
    - `R::Matrix`: The data to project, time x neurons
    - `B::Matrix`: The behavior data, time x properties
    # Returns
    - `Rproj::Matrix`: The projected behavior data, time x properties
    """
    function project_onto_behavior(em::AbstractMatrix, B::AbstractMatrix)
        lsq = llsq(Float64.(em), B);
        emp = (em * lsq[1:end-1,:]) .+ DIutils.arr.atleast2d(lsq[end, :])'
        lsq_coef = lsq[1:end-1,:]
        lsq_bias = lsq[end, :]
        return fit(PCA, emp, maxoutdim=3, pratio=1.0)
    end
    # """
    #         project_onto_behavior(R::Matrix, beh::DataFrame, beh_vars::Vector)
    # Project the data `R` onto the behavior variables `beh_vars` in `beh`.
    # # Arguments
    # - `R::Matrix`: The data to project, time x neurons
    # - `beh::DataFrame`: The behavior data, time x properties
    # # Returns
    # - `Rproj::Matrix`: The projected behavior data, time x properties
    # """
    # function project_onto_behavior(r::AbstractMatrix, beh::DataFrame, times, beh_vars::Vector)
    #     B = Matrix(beh[times, beh_vars])
    #     # Matrix(beh[times, beh_vars]) * r' * inv(r * r') * r
    #     r * B' * pinv(B * B') * B
    # end
    

    """
        register_two_manis(ref_emp, this_emp)
    Register two manifolds together using canonical correlation analysis.
    # Arguments
    - `ref_emp::EmbeddingProjection`: The reference embedding projection
    - `this_emp::EmbeddingProjection`: The embedding projection to register
    - `B::Matrix`: The behavior matrix, if you want to project onto behavior
    - `beh::DataFrame`: The behavior data, if you want to color plots by behavior
    - `filter_quantiles`: Whether to filter the data by quantiles
    # Returns
    - `M::CCA`: The CCA object
    - `Z1::Matrix`: The projection of the reference embedding onto the common space
    - `Z2::Matrix`: The projection of the embedding to register onto the common space
    - `C::Matrix`: The correlation between the projections
    """ 
    function pairwise_register(ref_emp::PCA, this_emp::PCA; B=nothing, beh=nothing,
            filter_quantiles=false, outdim=3)
        pairwise_register(ref_emp.proj, this_emp.proj; B=B, beh=beh,
            filter_quantiles=filter_quantiles, colorby=colorby, ploton=ploton, outdim=outdim)
    end
    function pairwise_register(ref_emp::Matrix, this_emp::Matrix; B=nothing, beh=nothing,
            filter_quantiles=false, outdim=3)
            # Cannonical corr
            XX = (ref_emp)'
            YY=  (this_emp)'
            if B !== nothing
                XX = project_onto_behavior(XX, B).proj
                YY = project_onto_behavior(YY, B).proj
            end
            if filter_quantiles
                xx = get_quantile_filtered(XX)
                yy = get_quantile_filtered(YY)
            else
                xx = XX
                yy = YY
            end
            # Randomly sample the same number of points
            n = min(size(xx,2), size(yy,2))
            xx = xx[:, sample(1:size(xx,2), n, replace=false)]
            yy = yy[:, sample(1:size(yy,2), n, replace=false)]
            M = fit(CCA, xx, yy, outdim=outdim)
            # Bring both into a common space
            a = xprojection(M);
            b = yprojection(M);
            Z1 = StatsAPI.predict(M, (XX), :x)
            Z2 = StatsAPI.predict(M, (YY), :y)
            return (;M, Z1, Z2)
    end


end
