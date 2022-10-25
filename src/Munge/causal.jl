module causal
    using CausalityTools
    using Entropies

    export get_est_preset
    function get_est_preset(esttype::Symbol)
        if esttype == :binned
            params = (;bins = 9, horizon=1:3000)
            est = VisitationFrequency(RectangularBinning(params.bins))
            @info "Starting binned estimator, "
        elseif esttype == :symbolic
            params = (m = 15, τ = 1)
            #params = (m = 100, τ = 4)
            est = SymbolicPermutation(;params...)
            params = (;params..., horizon=1:3000)
        else
            @warn "not recog"
        end
        (;est, params)
    end
    
    export global_predictive_asymmetry
    function global_predictive_asymmetry(embeddingX::DataSet,
            embeddingY::DataSet, est; params...)
        Thrads.@spawn CausalityTools.predictive_asymmetry(embeddingX, 
                                            embeddingY, 
                                            est,
                                            params.horizon)
    end
    function global_predictive_asymmetry(embeddingX::AbstractArray,
                embeddingY::AbstractArray, est; params)
        global_predictive_asymmetry(
             DataSet(embeddingX), DataSet(embeddingY), est; params...)
    end
    function global_predictive_asymmetry(embeddingX::AbstractDict,
                                        embeddingY::AbstractDict, est;
                                        results=Dict(),
                                        params...)
        keysX, keysY = keys(embeddingX), keys(embeddingY) 
        keysX = [(k,v) for (k,v) in keysX if k!=pairdim]
        keysY = [(k,v) for (k,v) in keysX if k!=pairdim]
        kz = intersect(keysX, keysY)

        for k in kz
            # For right now, these pairing dimensions are hard-coded
            k1=bestpartialmatch(keys(embeddingX), (;k..., area=:ca1))
            k2=bestpartialmatch(keys(embeddingX), (;k..., area=:pfc))
            results[k] = global_predictive_asymmetry(embeddingX[k1],
                                                     embeddingY[k2],
                                                     est;params...)
        end
    end

    function local_predictive_asymmetry()
    end
end
