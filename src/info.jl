module info
    using Entropies: Probabilities
    skipnan(x) = Iterators.filter(!isnan, x)

    """
    `spatial_information`

    computes the `spatial_information` 

    see yartsev dotson 2021 supp
    """
    function spatial_information(fields::Dict, 
            behProb::Union{AbstractArray, Probabilities})
        if behfield isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behfield))))
        end
        I = copy(fields)
        for i in keys(fields)
            I[i] = spatial_information(fields)
        end
        return I
    end
    function spatial_information(fields::AbstractArray, 
            behProb::Union{AbstractArray, Probabilities})
        if behfield isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behfield))))
        end
        I = copy(fields)
        fields = collect(skipnan(vec(fields)))
        R = fields ./ mean(fields)
        I = sum( behProb .* R .* log2.(R) )
        return I
    end

end
export info
