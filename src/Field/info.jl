module info
    using Entropies: Probabilities
    using NaNStatistics
    import ..Field
    using Infiltrator

    skipnan(x) = Iterators.filter(!isnan, x)
    function _convert_to_prob(behProb::AbstractArray)
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
    end

    """
    `information`

    computes the `information` of a receptive field ... general dispatcher for
    a number of information metrics

    see yartsev dotson 2021 supp
    """
    function information(F::Field.ReceptiveField; method=:spatialinformation)
        return information(F.rate, F.occ.prob; method)
    end
    function information(F::NamedTuple; method=:spatialinformation)
        return information(F.Rₕ, F.occR; method)
    end
    function information(F::Dict, behProb::AbstractArray, method=:spatialinformation)
        behProb = info._convert_to_prob(behProb)
        information(F, behProb; method)
    end
    function information(F::Dict, behProb::Probabilities; method=:spatialinformation)
        I = Dict()
        method = eval(method)
        for i in keys(F)
            I[i] = method(F[i], behProb)
        end
        return I
    end
    function information(F::Array, behProb::Probabilities; method=:spatialinformation)
        eval(method)(F, behProb)
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
        
        beh, spikes = Load.keep_overlapping_times(beh, spikes)


        

    end

end
export info
