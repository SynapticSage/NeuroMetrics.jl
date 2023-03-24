module causal

    using Infiltrator

    using ProgressMeter
    using CausalityTools
    using Entropies
    using CausalityTools: Dataset
    using Infiltrator
    using DataFrames
    using ShiftedArrays
    using DataStructures: OrderedDict
    using StaticArrays
    import DelayEmbeddings
    using ThreadSafeDicts
    using JLD2
    using SoftGlobalScope
    using DrWatson
    using ArgParse

    import ..Munge
    import ..Munge.manifold: make_embedding_df, EmbeddingFrameFetch
    import DIutils
    import DIutils: Table
    import DIutils.namedtup: bestpartialmatch

    function argparse(args=nothing;return_parser=false)
        parser = Munge.manifold.parse(return_parser=true)
        if return_parser
            parser
        else
            args === nothing ? parse_args(parser) : parse_args(parser, args)
        end
    end

    export get_est_preset
    function get_est_preset(esttype::Symbol; 
            horizon=1:300, bins = [0.25,0.25,0.25],
            m=15, τ=1, other_params...
        )
        if esttype == :binned
            params = (;bins, horizon)
            est = VisitationFrequency(RectangularBinning(params.bins))
            @info "Starting binned estimator, "
        elseif esttype == :symbolic
            #params = (m = 100, τ = 4)
            params = (m, τ)
            est = SymbolicPermutation(;params...)
            params = (;params..., horizon)
        else
            @error "not recog"
        end
        params = (;params..., other_params...)
        (;est, params)
    end

    export get_paired_emeddings
    function get_paired_embeddings(embeddings::AbstractDict,
                                  valuesX, valuesY; prop=:dataset)::Dict

        kz = keys(embeddings)
        kz = unique(DIutils.namedtup.pop.(kz,[:dataset]))

        results = Dict{NamedTuple, Tuple}()
        for k in kz
            # For right now, these pairing dimensions are hard-coded
            kX, kY = (;k..., dataset=valuesX), 
                     (;k..., dataset=valuesY)
            kX = bestpartialmatch(keys(embeddings),kX)
            kY = bestpartialmatch(keys(embeddings),kY)
            if valuesX != valuesY
                @assert kX != kY "kX and kY should be different if valuesX != valuesY"
            end
            results[k] = (embeddings[kX], embeddings[kY])
        end
        results
    end

    
    ## ---------- GLOBAL WITH AN ESTIMATOR ------------------
    export predictiveasymmetry
    function predictiveasymmetry(embeddingX::Dataset,
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
    function predictiveasymmetry(embeddingX::AbstractArray,
                embeddingY::AbstractArray, est; params...)
        if size(embeddingX,2) > size(embeddingX,1)
            @assert size(embeddingY,2) > size(embeddingY,1)
            embeddingX, embeddingY = embeddingX', embeddingY'
        end
        causal.predictiveasymmetry(
             Dataset(embeddingX), Dataset(embeddingY), est; params...)
    end
    function predictiveasymmetry(pairedembeddings::Dict, 
            est::ProbabilitiesEstimator; 
            params...)
        Dict(k=> predictiveasymmetry(v[1],v[2],est;params...)
                 for (k,v) in pairedembeddings)
    end
    function predictiveasymmetry!(checkpoint::AbstractDict, 
            pairedembeddings::Dict, est; params...)
        for (k,v) in pairedembeddings
            if k ∉ keys(checkpoint)
                push!(checkpoint, k=>predictiveasymmetry(v[1],v[2],est;params...))
            end
        end
    end

    ## ---------- GLOBAL WITHOUT AN ESTIMATOR ------------------
    function predictiveasymmetry(embeddingX::Dataset,
            embeddingY::Dataset; thread=true, params...)

        # Subset by a condition?
        if :condition_inds ∈ keys(params) && 
            :s ∈ keys(params) && 
            :inds_of_t ∈ keys(params)
            @info "indexing" params[:s]
            subset = params[:condition_inds][collect(params[:inds_of_t][params[:s]])]
            embeddingX, embeddingY = embeddingX[subset], embeddingY[subset]
        else
            #@info "not indexing"
        end
        
        # Embed into univariate space
        if !isempty(embeddingX)
            est  = RectangularBinning(params[:binning])
            uniX = encode_dataset_univar(embeddingX, est)
            uniY = encode_dataset_univar(embeddingY, est)
            M = Int64(max(maximum(uniX),maximum(uniY)))
            est = VisitationFrequency(RectangularBinning(M))
            #@info est
        end

        if !isempty(embeddingX) && maximum(params[:horizon]) > length(uniX)
            @error "Horizon too long"
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
            CausalityTools.predictive_asymmetry(uniX, 
                                                uniY, 
                                                est,
                                                params[:horizon])
        end
    end
    function predictiveasymmetry(embeddingX::AbstractArray,
                embeddingY::AbstractArray; params...)
        if size(embeddingX,2) > size(embeddingX,1)
            @assert size(embeddingY,2) > size(embeddingY,1)
            embeddingX, embeddingY = embeddingX', embeddingY'
        end
        predictiveasymmetry(
             Dataset(embeddingX), Dataset(embeddingY); params...)
    end
    function predictiveasymmetry(pairedembeddings::Dict; 
            params...)
        Dict(k=> predictiveasymmetry(v[1],v[2];k...,params...)
                 for (k,v) in pairedembeddings)
    end
    export predictiveasymmetry!
    function predictiveasymmetry!(checkpoint::AbstractDict,
            pairedembeddings::Dict; params...)
        for (k,v) in pairedembeddings
            if k ∉ keys(checkpoint)
                push!(checkpoint, k=>predictiveasymmetry(v[1],v[2];k...,params...))
            elseif k ∈ keys(checkpoint) # if we re-encounter a key! (possible checkpoint computation)
                previous = checkpoint[k]
                # We are reencountering a key, is it a task, if so fetch it
                if previous isa Task
                    checkpoint[k] = previous = fetch(checkpoint[k])
                end
                # If we have a vector, is it the expected length for a given horizon; if so
                # skip it, else compute the missing horizon
                if previous isa Vector
                    if length(previous) == length(params[:horizon])
                        @info "predictiveasymmetry reencoutered key, but already computed, skipping" 
                        continue
                    else
                        @info "predictiveasymmetry reencountered key, and extending horizon" length(previous) length(params[:horizon])
                        newhorizon = setdiff(collect(params[:horizon]), 1:length(params[:horizon]))
                        newparams = (;params..., horizon=newhorizon)
                        push!(checkpoint, k=>predictiveasymmetry(v[1],v[2];
                                                                 k...,newparams...))
                    end
                else
                    @warn "predictiveasymmetry may have reencountered a key which prevous failed" key=k
                end
            end
        end
    end

    ## ---------- CONDITIONAL WITHOUT AN ESTIMATOR ------------------
    export conditional_pred_asym
    function conditional_pred_asym(checkpoint::AbstractDict, 
            pairedembeddings::Dict, data::DataFrame, data_vars; 
            groups=nothing, params...)

        @info "conditional_pred_asym" data_vars groups 

        # Figure out the possible permutations of the data-variables user wants to
        # condition on
        data_vars = replace(hcat([data[!,var] for var in data_vars]...),
                            NaN=>-1,missing=>-1)
        groupinds = DIutils.findgroups(data_vars)
        groupsfulllist    = OrderedDict(data_vars[findfirst(groupinds.==k),:]=>k
                          for k in unique(groupinds))

        # Enumerate list, not given by user
        if groups === nothing
            groups = collect(keys(groupsfulllist))
        end

        for group in groups
            @info group
            # If key not in list, then  skip it
            if group ∉ keys(groupsfulllist)
                @info "runconditional: group=$group ∉ keys(groupsfulllist)"
                continue
            end
            # Figure out which indices are our group of interest
            condition_inds = groupsfulllist[group] .== groupinds
            # If key doesn't exist, make a new dict there
            if group ∉ keys(checkpoint)
                checkpoint[group] = Dict()
            end
            # If dict type incorrect, fix it
            if !(Task <: valtype(checkpoint[group])) || 
                !(Task <: valtype(checkpoint[group]))
                kt, vt = keytype(checkpoint), Union{Task,Vector,valtype(checkpoint[group])}
                checkpoint[group] = Dict{kt}{vt}(checkpoint[group])
            end
            predictiveasymmetry!(checkpoint[group], pairedembeddings; condition_inds, params...)
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
        convert(Vector{Int16}, DIutils.searchsortednearest.([sout], vals))
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

    export obtain_triggered_causality_binavg
    """
        obtain_triggered_causality(em, beh, props, params; grid_kws=nothing, 
            savefile,
            grd=(println("default binning compute");
                 binning.get_grid(beh, props; grid_kws...)),
            checkpoint=ThreadSafeDict{NamedTuple, Bool}())

    Gets the triggered causality, triggered by samples entering a bin. 
    """
    function obtain_triggered_causality_binavg(em::DataFrame, beh::DataFrame, 
            props::Vector, params::NamedTuple; 
            grid_kws=nothing, savefile,
            grd=(println("default binning compute");
                 binning.get_grid(beh, props; grid_kws...)),
            checkpoint=ThreadSafeDict{NamedTuple, Bool}(),
            checkpointmod = 1, 
            inttype::Type=Int32, floattype::Type=Float32, nan=NaN32,
        )

        @info "obtain_local_binned_measure()" props savefile params params.window
        ## ------------------------------
        ## Setup frames and zone triggers
        ## ------------------------------
        #propstime = ["time",props...]
        EFF = EmbeddingFrameFetch(em, :area, beh, props; 
                                  ordering=Dict(:area=>[:ca1,:pfc]))

        ## ------------------------
        ## Setup stores for results
        ## ------------------------
        ca1pfc = fill(nan,size(grd)..., length(EFF), params[:horizon].stop)
        pfcca1 = fill(nan,size(grd)..., length(EFF), params[:horizon].stop)
        ca1pfcσ² = fill(nan,size(grd)..., length(EFF), params[:horizon].stop)
        pfcca1σ² = fill(nan,size(grd)..., length(EFF), params[:horizon].stop)
        counts   = fill(inttype(0),size(grd)..., length(EFF))
        nancounts = fill(inttype(0), size(grd)..., length(EFF))

        if isfile(savefile)
            @info "obtain_local_binned_measure, found save file, loading" savefile
            storage    = JLD2.jldopen(savefile, "r")
            checkpoint = merge(storage["checkpoint"], checkpoint)
            ca1pfc = "ca1pfc" ∈ keys(storage) ? storage["ca1pfc"]  : ca1pfc
            pfcca1 = "pfcca1" ∈ keys(storage) ? storage["pfcca1"] : pfcca1
            ca1pfcσ² = "ca1pfcσ²" ∈ keys(storage) ? storage["ca1pfcσ²"] : ca1pfcσ²
            pfcca1σ² = "pfcca1σ²" ∈ keys(storage) ? storage["pfcca1σ²"] : pfcca1σ²
            counts = "counts" in keys(storage) ? storage["counts"] : counts
            nancounts = "counts" in keys(storage) ? storage["nancounts"] : nancounts
            close(storage)
        end
        

        ## ------------------------
        ## Compute
        ## ------------------------
        #(e,eff) = first(enumerate(EFF))

        @info "How much ground to cover" length(EFF) length(first(EFF))
        @showprogress "frames" for (iEmbed,eff) in enumerate(EFF[1:length(EFF)])
            @info iEmbed 
            zones = Munge.triggering.get_triggergen(params.window, grd, eff...)
            #(idx,data) = zones[50]
            for (idx,data) in zones
                if isempty(data); continue; end
                checkpointkey = (;idx, iEmbed)
                if checkpointkey ∈ keys(checkpoint) && checkpoint[checkpointkey]
                    continue
                else
                    checkpoint[checkpointkey] = false
                end
                count = inttype(0)
                ncount = inttype(0)
                cpv = view(ca1pfc, idx..., iEmbed, :)
                pcv = view(pfcca1, idx..., iEmbed, :)
                cpvσ² = view(ca1pfcσ², idx..., iEmbed, :)
                pcvσ² = view(pfcca1σ², idx..., iEmbed, :)
                Threads.@threads for D in data
                    if isempty(D); continue; end
                    ca1, pfc = map(d->transpose(hcat(d.data...)), D)
                    try
                        cpt = causal.predictiveasymmetry(ca1, pfc; params...)
                        pct = causal.predictiveasymmetry(pfc, ca1; params...)
                        if any(isnan.(cpt)) || any(isnan.(pct))
                            ncount .+= 1
                            continue
                        end
                        cpt,pct = convert(Vector{floattype}, cpt),
                                  convert(Vector{floattype}, pct)
                        #if all(isnan.(cpv)) && all(isnan.(pcv))
                        #    cpv .= 0
                        #    pcv .= 0
                        #end
                        cpv .+= cpt
                        pcv .+= pct
                        cpvσ² .+= cpt.^2
                        pcvσ² .+= pct.^2
                        
                    catch err
                        showerror(stdout, err)
                    end
                    count += 1
                end
                ca1pfc[idx..., iEmbed, :]   ./= count
                pfcca1[idx..., iEmbed, :]   ./= count
                pfcca1σ²[idx..., iEmbed, :] ./= count
                ca1pfcσ²[idx..., iEmbed, :] ./= count
                counts[idx...,   iEmbed] = count
                nancounts[idx..., iEmbed] = ncount
                checkpoint[checkpointkey] = true
            end
            if mod(iEmbed,checkpointmod) == 0
                storage    = jldopen(savefile, "a")
                #@info "showing store" storage
                ["ca1pfc","pfcca1","checkpoint","counts", "ca1pfcσ²","pfcca1σ²"]
                "ca1pfc" in keys(storage) ? delete!(storage,"ca1pfc") : nothing
                storage["ca1pfc"] = ca1pfc
                "pfcca1" in keys(storage) ? delete!(storage,"pfcca1") : nothing
                storage["pfcca1"] = pfcca1
                "counts" in keys(storage) ? delete!(storage,"counts") : nothing
                storage["counts"] = counts
                "checkpoint" in keys(storage) ? delete!(storage,"checkpoint") : nothing
                storage["checkpoint"] = checkpoint
                "ca1pfcσ²" in keys(storage) ? delete!(storage,"ca1pfcσ²") : nothing
                storage["ca1pfcσ²"] = ca1pfcσ²
                "pfcca1σ²" in keys(storage) ? delete!(storage,"pfcca1σ²") : nothing
                storage["pfcca1σ²"] = pfcca1σ²
                close(storage)
            end
        end
        savefile
    end

    # ANALYSIS SAVE FILES
    """
        find_horizon_values

    obtains all horizon values that are present in files in the munge/causal savefolder
    optionally contraining to values less than or equal to `constraint`
    """
    function _find_horizon_values(;constraint::Union{Int,Nothing}=nothing)
        files = readdir(datadir("manifold","causal")) 
        search=",horizon=>"
        # Cut out the horizon data from the string names
        inds = findfirst.([search], files)
        data = [file[last(ind)+1:end] for (file,ind) in zip(files,inds)
               if ind !== nothing]
        inds = findfirst.(",", data)
        inds = [(ind === nothing ? findfirst("_",datum) : ind) for (datum,ind) in zip(data,inds)]
        data = [file[begin:first(ind)-1] for (file,ind) in zip(data,inds)]
        data = replace.(data, "-"=>":")
        data = unique(last.([Base.eval(instr) for instr in Meta.parse.(data)]))
        constraint === nothing ? data : data[data .<= constraint]
    end
    export get_alltimes_savefile
    """
        get_trigger_savefile

    Obtains the savefile name for given `params`
    """
    function get_alltimes_savefile(animal, day, N; params=(;), allowlowerhorizon=false) 
        params = pop!(params,:thread)
        tagstr = "$animal$day.$(N)seg"
        @info "param order " keys(params)
        filename = ""
        # allowhorizon means we search for multiple possible filenames, as long as the other params are the same and horizon less than what user asked for
        if allowlowerhorizon
            @info "Looking for file with horizon <= $(last(params[:horizon]))"
            params = OrderedDict(pairs(params))
            horizons = _find_horizon_values(;constraint=last(params[:horizon]))
            possparams = [(setindex!(copy(params), 1:horizon, :horizon))
                      for horizon in horizons]
            paramstrs = DIutils.namedtup.tostring.(NamedTuple.(possparams),keysort=true)
            files = [datadir("manifold", "causal", "pa_cause_$(paramstr)_$tagstr.jld2")
                    for paramstr in paramstrs]
            exists = isfile.(files)
            if any(exists)
                horizons, files = horizons[exists], files[exists]
                filename = files[argmax(horizons)]
            else
                filename = ""
            end
        end
        if isempty(filename)
            paramstr = DIutils.namedtup.tostring(NamedTuple(params); keysort=true)
            filename = datadir("manifold","causal", "pa_cause_$(paramstr)_$tagstr.jld2")
        end
        filename
    end

    export load_alltimes_savefile
    function load_alltimes_savefile(animal, day, N; params=(;))
        jldopen(get_alltimes_savefile(animal,day,N;params), "r")
    end
    
    export get_trigger_savefile
    """
        get_trigger_savefile

    Obtains the savefile name for given `props` and `params`
    """
    function get_trigger_savefile(animal, day, N, props; params=(;)) 
        tagstr = "$animal$day.$(N)seg"
        paramstr = DIutils.namedtup.tostring(pop!(params,:thread))
        datadir("manifold","causal",
                "local_grid_cause_props=$(join(props,","))_$(paramstr)_$tagstr.jld2")
    end
    """
        load_trigger_savefile

    Handles the loading of a trigger save file. See docfile for `get_trigger_savefile`
    """
    function load_trigger_savefile(animal, day, N, props; params=(;))
        savefile = get_trigger_savefile(animal, day, N, props; params)
        @assert isfile(savefile) "file is missing"
        jldopen(savefile,"r")
    end

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
