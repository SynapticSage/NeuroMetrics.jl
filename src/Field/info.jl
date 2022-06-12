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
        I = nansum( behProb .* R .* log2.(R) )
        return I
    end

end
export info
