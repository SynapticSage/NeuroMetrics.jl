module causal

    import Utils

    using CausalityTools
    using Entropies
    using CausalityTools: Dataset
    using Infiltrator

    export get_est_preset
    function get_est_preset(esttype::Symbol)
        if esttype == :binned
            params = (;bins = [0.25,0.25,0.25], horizon=1:2500)
            est = VisitationFrequency(RectangularBinning(params.bins))
            @info "Starting binned estimator, "
        elseif esttype == :symbolic
            params = (m = 15, τ = 1)
            #params = (m = 100, τ = 4)
            est = SymbolicPermutation(;params...)
            params = (;params..., horizon=1:2500)
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

    
    export global_predictive_asymmetry
    function global_predictive_asymmetry(embeddingX::Dataset,
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
    function global_predictive_asymmetry(embeddingX::AbstractArray,
                embeddingY::AbstractArray, est; params...)
        global_predictive_asymmetry(
             Dataset(embeddingX), Dataset(embeddingY), est; params...)
    end
    function global_predictive_asymmetry(pairedembeddings::Dict, est; 
            params...)
        Dict(
             k=> global_predictive_asymmetry(v[1],v[2],est;params...)
                 for (k,v) in pairedembeddings)
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
