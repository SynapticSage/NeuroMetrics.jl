module causal

    using Infiltrator
    import Utils

    using CausalityTools
    using Entropies
    using CausalityTools: Dataset
    using Infiltrator
    using DataFrames
    using ShiftedArrays
    using DataStructures: OrderedDict
    using StaticArrays
    import Table
    import DelayEmbeddings
    #using PyCall
    #@pyimport PyIF

    import Munge.manifold: make_embedding_df, EmbeddingFrameFetch

    export get_est_preset
    function get_est_preset(esttype::Symbol, horizon=1500)
        if esttype == :binned
            params = (;bins = [0.25,0.25,0.25], horizon=1:horizon)
            est = VisitationFrequency(RectangularBinning(params.bins))
            @info "Starting binned estimator, "
        elseif esttype == :symbolic
            params = (m = 15, τ = 1)
            #params = (m = 100, τ = 4)
            est = SymbolicPermutation(;params...)
            params = (;params..., horizon=1:horizon)
        else
            @warn "not recog"
        end
        (;est, params)
    end

    export get_paired_emeddings
    function get_paired_embeddings(embeddings::AbstractDict,
                                  valuesX, valuesY; prop=:area)::Dict

        kz = keys(embeddings)
        kz = unique(Utils.namedtup.pop.(kz,[:area]))

        results = Dict{NamedTuple, Tuple}()
        for k in kz
            # For right now, these pairing dimensions are hard-coded
            kX, kY = (;k..., area=valuesX), (;k..., area=valuesY)
            kX = Utils.bestpartialmatch(keys(embeddings),kX)
            kY = Utils.bestpartialmatch(keys(embeddings),kY)
            results[k] = (embeddings[kX], embeddings[kY])
        end
        results
    end

    
    ## ---------- GLOBAL WITH AN ESTIMATOR ------------------
    export predictive_asymmetry
    function predictive_asymmetry(embeddingX::Dataset,
            embeddingY::Dataset, est; thread=true, params...)
        if thread
            Threads.@spawn CausalityTools.predictive_asymmetry(
                                                embeddingX, 
                                                embeddingY, 
                                                est,
                                                params[:horizon])
        else
            CausalityTools.predictive_asymmetry(embeddingX, 
                                                embeddingY, 
                                                est,
                                                params[:horizon])
        end
    end
    function predictive_asymmetry(embeddingX::AbstractArray,
                embeddingY::AbstractArray, est; params...)
        predictive_asymmetry(
             Dataset(embeddingX), Dataset(embeddingY), est; params...)
    end
    function predictive_asymmetry(pairedembeddings::Dict, 
            est::ProbabilitiesEstimator; 
            params...)
        Dict(k=> predictive_asymmetry(v[1],v[2],est;params...)
                 for (k,v) in pairedembeddings)
    end
    function predictive_asymmetry!(checkpoint::AbstractDict, 
            pairedembeddings::Dict, est; params...)
        for (k,v) in pairedembeddings
            if k ∉ keys(checkpoint)
                push!(checkpoint, k=>predictive_asymmetry(v[1],v[2],est;params...))
            end
        end
    end

    ## ---------- GLOBAL WITHOUT AN ESTIMATOR ------------------
    function predictive_asymmetry(embeddingX::Dataset,
            embeddingY::Dataset; thread=true, params...)

        # Subset by a condition?
        if :condition_inds ∈ keys(params) && 
            :s ∈ keys(params) && 
            :inds_of_t ∈ keys(params)
            @info "indexing" params[:s]
            subset = params[:condition_inds][collect(params[:inds_of_t][params[:s]])]
            embeddingX, embeddingY = embeddingX[subset], embeddingY[subset]
        else
            @info "not indexing"
        end
        
        # Embed into univariate space
        if !isempty(embeddingX)
            est  = RectangularBinning(params[:binning])
            uniX = encode_dataset_univar(embeddingX, est)
            uniY = encode_dataset_univar(embeddingY, est)
            M = Int64(max(maximum(uniX),maximum(uniY)))
            est = VisitationFrequency(RectangularBinning(M))
            @info est
        end

        if isempty(embeddingX)
            missing
        elseif thread
            Threads.@spawn @time CausalityTools.predictive_asymmetry(
                                                uniX, 
                                                uniY, 
                                                est,
                                                params[:horizon])
        else
            @time CausalityTools.predictive_asymmetry(uniX, 
                                                uniY, 
                                                est,
                                                params[:horizon])
        end
    end
    function predictive_asymmetry(embeddingX::AbstractArray,
                embeddingY::AbstractArray; params...)
        predictive_asymmetry(
             Dataset(embeddingX), Dataset(embeddingY); params...)
    end
    function predictive_asymmetry(pairedembeddings::Dict; 
            params...)
        Dict(k=> predictive_asymmetry(v[1],v[2];k...,params...)
                 for (k,v) in pairedembeddings)
    end
    export predictive_asymmetry!
    function predictive_asymmetry!(checkpoint::AbstractDict,
            pairedembeddings::Dict; params...)
        for (k,v) in pairedembeddings
            if k ∉ keys(checkpoint)
                push!(checkpoint, k=>predictive_asymmetry(v[1],v[2];k...,params...))
            end
        end
    end

    ## ---------- CONDITIONAL WITHOUT AN ESTIMATOR ------------------
    export conditional_pred_asym
    function conditional_pred_asym(checkpoint::AbstractDict, 
            pairedembeddings::Dict, data::DataFrame, data_vars=[:cuemem, :correct]; 
            groups=nothing, params...)

        data_vars = replace(hcat([data[!,var] for var in data_vars]...),NaN=>-1,missing=>-1)
        groupinds = Utils.findgroups(data_vars)
        groupsfulllist    = OrderedDict(data_vars[findfirst(groupinds.==k),:]=>k
                          for k in unique(groupinds))
        if groups === nothing
            groups = groupsfulllist
        end
        #groups    = [[1,1],[0,1],[1,0],[0,0]]

        for group in groups
            @info group
            condition_inds = groupsfulllist[group] .== groupinds
            if group ∉ keys(checkpoint)
                checkpoint[group] = Dict()
            end
            predictive_asymmetry!(checkpoint[group], pairedembeddings; condition_inds, params...)
        end

    end

    ##-------------------------------------------------------------
    ## STOLEN FROM ENTROPIES.JL  --- Missing these off develop ---
    ##-------------------------------------------------------------
    """
    Encoding
    The supertype for all encoding schemes, i.e. ways of encoding input data onto some
    discrete set of possibilities/[outcomes](@ref). Implemented encoding schemes are
    - [`OrdinalPatternEncoding`](@ref).
    - [`GaussianCDFEncoding`](@ref).
    - [`RectangularBinEncoding`](@ref).
    Used internally by the various [`ProbabilitiesEstimator`](@ref)s to map input data onto
    outcome spaces, over which probabilities are computed.
    """
    abstract type Encoding end
    """
        outcomes(x, scheme::Encoding) → Vector{Int}
        outcomes!(s, x, scheme::Encoding) → Vector{Int}
    Map each`xᵢ ∈ x` to a distinct outcome according to the encoding `scheme`.
    Optionally, write outcomes into the pre-allocated symbol vector `s` if the `scheme`
    allows for it. For usage examples, see individual encoding scheme docstrings.
    See also: [`RectangularBinEncoding`](@ref), [`GaussianCDFEncoding`](@ref),
    [`OrdinalPatternEncoding`](@ref).
    """
    function outcomes(x::X, ::Encoding) where X
        throw(ArgumentError("`outcomes` not defined for input data of type $(X)."))
    end
    """
    RectangularBinEncoding <: Encoding
    RectangularBinEncoding(x, binning::RectangularBinning)
    RectangularBinEncoding(x, binning::FixedRectangularBinning)
    Find the minima along each dimension, and compute appropriate
    edge lengths for each dimension of `x` given a rectangular binning.
    Put them in an `RectangularBinEncoding` that can be then used to map points into bins
    via [`outcomes`](@ref).
    See also: [`RectangularBinning`](@ref), [`FixedRectangularBinning`](@ref).
    """
    struct RectangularBinEncoding{B, M, E} <: Encoding
        binning::B # either RectangularBinning or FixedRectangularBinning
        mini::M # fields are either static vectors or numbers
        edgelengths::E
    end
    function encode_as_bin(point, b::RectangularBinEncoding)
        (; mini, edgelengths) = b
        # Map a data point to its bin edge
        return floor.(Int, (point .- mini) ./ edgelengths)
    end
    function outcomes(x::Dataset, b::RectangularBinEncoding)
        return map(point -> encode_as_bin(point, b), x)
    end
    function RectangularBinEncoding(x::Dataset{D,T}, b::RectangularBinning;
        n_eps = 2) where {D, T}
        # This function always returns static vectors and is type stable
        ϵ = b.ϵ
        mini, maxi = DelayEmbeddings.minmaxima(x)
        v = ones(SVector{D,T})
        if ϵ isa Float64 || ϵ isa AbstractVector{<:AbstractFloat}
            edgelengths = ϵ .* v
        elseif ϵ isa Int || ϵ isa Vector{Int}
            edgeslengths_nonadjusted = @. (maxi - mini)/ϵ
            # Just taking nextfloat once here isn't enough for bins to cover data when using
            # `encode_as_bin` later, because subtraction and division leads to loss
            # of precision. We need a slightly bigger number, so apply nextfloat twice.
            edgelengths = nextfloat.(edgeslengths_nonadjusted, n_eps)
        else
            error("Invalid ϵ for binning of a dataset")
        end
        RectangularBinEncoding(b, mini, edgelengths)
    end
    # ----------- END OF STEAL --------------------------
    
    export encode_dataset
    function encode_dataset(x::Dataset, ϵ::RectangularBinning)
        enc =  RectangularBinEncoding(x, ϵ)
        outcomes(x, enc)
    end
    export encode_dataset_univar
    function encode_dataset_univar(x::Dataset, ϵ::RectangularBinning)
        out = encode_dataset(x, ϵ)
        out = Matrix(hcat([collect(x) for x in out]...)')
        M = maximum(out)
        vals = vec(sum(Matrix(out .* (10 .^ (1:size(out,2)))'),dims=2))
        sout = sort(unique(vals))
        convert(Vector{Int16}, Utils.searchsortednearest.([sout], vals))
    end



    export predictive_asymmetry_pyif
    function predictive_asymmetry_pyif(x, y, horizon, pyif::Module; 
            embedding=1, knearest=1, gpu=false)

        PA,AP=[],[]
        for h in horizon
            yy = replace(ShiftedArray(y, (h)), missing=>0)
            push!(PA,
                  pyif.te_compute.te_compute(x,y,
                    embedding=embedding,k=knearest,GPU=gpu)
                 )
            push!(AP,
                  pyif.te_compute.te_compute(y,x,
                        embedding=embedding,k=knearest,GPU=gpu)
                 )
        end
        cumsum(PA) - cumsum(AP)
    end

    #function global_predictive_asymmetry(embeddingX::AbstractDict,
    #                                    embeddingY::AbstractDict, est;
    #                                    results=Dict(),
    #                                    params...)
    #    keysX, keysY = keys(embeddingX), keys(embeddingY) 
    #    keysX = [(k,v) for (k,v) in keysX if k!=pairdim]
    #    keysY = [(k,v) for (k,v) in keysX if k!=pairdim]
    #    kz = intersect(keysX, keysY)

    #    for k in kz
    #        # For right now, these pairing dimensions are hard-coded
    #        k1=bestpartialmatch(keys(embeddingX), (;k..., area=:ca1))
    #        k2=bestpartialmatch(keys(embeddingX), (;k..., area=:pfc))
    #        results[k] = global_predictive_asymmetry(embeddingX[k1],
    #                                                 embeddingY[k2],
    #                                                 est;params...)
    #    end
    #end

    #function local_predictive_asymmetry()
    #end
end
