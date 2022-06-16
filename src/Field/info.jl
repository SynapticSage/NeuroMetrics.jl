module info
    using Entropies: Probabilities
    using NaNStatistics
    skipnan(x) = Iterators.filter(!isnan, x)

    """
    `information`

    computes the `information` of a receptive field

    see yartsev dotson 2021 supp
    """
    function information(F::NamedTuple)
        return info.information(F.Rₕ, F.occR)
    end
    function information(F::Dict, 
            behProb::Union{AbstractArray, Probabilities})
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
        I = Dict()
        for i in keys(F)
            I[i] = info.information(F[i], behProb)
        end
        return I
    end
    function information(field::AbstractArray, 
            behProb::Union{AbstractArray, Probabilities})
        if behProb isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behProb))))
        end
        #I = copy(field)
        field = collect(skipnan(vec(field)))
        R = field ./ nanmean(field)
        # paper r = ∑ pᵢ * rᵢ
        # what I've written here: r = ∑ R_occᵢ / N = μ(R_occᵢ), not the actual rate, but occ norm rate
        I = nansum( behProb .* R .* log2.(R) )
        # what i've written here: ∑ p(state) * r(state)/r̅ * log₂( r(state) / r̅ )
        # paper : ∑ pᵢ * ( rᵢ/r̅ ) * log₂( rᵢ / r̅ )
        return I
    end

    function mutualinformation()
    end

end
export info
