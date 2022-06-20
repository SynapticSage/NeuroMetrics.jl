module info
    using Entropies: Probabilities
    using NaNStatistics

    skipnan(x) = Iterators.filter(!isnan, x)
    function _convert_to_prop(behProb::AbstractArray)
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
    end

    """
    `information`

    computes the `information` of a receptive field

    see yartsev dotson 2021 supp
    """
    function information(F::NamedTuple; kws...)
        return info.information(F.Rₕ, F.occR; kws...)
    end
    function information(F::Dict, behProb::AbstractArray, method=:spatialinformation)
        behProb = _convert_to_prob(behProb)
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
    function bitsperspike(field::AbstractArray, behProb::Probabilities)
        field = collect(skipnan(vec(field)))
        R = field ./ nanmean(field)
        # paper r = ∑ pᵢ * rᵢ
        # what I've written here: r = ∑ R_occᵢ / N = μ(R_occᵢ), not the actual rate, but occ norm rate
        I = nansum( behProb .* R .* log2.(R) )
        # what i've written here: ∑ p(state) * r(state)/r̅ * log₂( r(state) / r̅ )
        # paper : ∑ pᵢ * ( rᵢ/r̅ ) * log₂( rᵢ / r̅ )
        return I
    end
    spatialinformation = bitsperspike
    function bitspersecond(field::AbstractArray, behProb::Probabilities)
        field = collect(skipnan(vec(field)))
        R = field ./ nanmean(field)
        # paper r = ∑ pᵢ * rᵢ
        # what I've written here: r = ∑ R_occᵢ / N = μ(R_occᵢ), not the actual rate, but occ norm rate
        I = nansum( behProb .* field .* log2.(R) )
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
    function mutualinformation(field::AbstractVector{<:Int}, behCount::AbstractArray{<:Int})
    end

end
export info
