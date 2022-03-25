module info
    using Entropies: Probabilities
    skipnan(x) = Iterators.filter(!isnan, x)

    """
    `spatial_information`

    computes the `spatial_information` 

    see yartsev dotson 2021 supp
    """
    function spatial_information(fields::Union{Dict,AbstractArray}, 
            behProb::Union{AbstractArray, Probabilities})
        if behfield isa AbstractArray
            behProb = Probabilities(collect(skipnan(vec(behfield))))
        end
        I = copy(fields)
        if fields isa Dict
            for i in keys(fields)
                I[i] = spatial_information(fields)
            end
        else
            fields = collect(skipnan(vec(fields)))
            R = fields ./ mean(fields)
            I = sum( behProb .* R .* log2.(R) )
        end
        return I
    end

    function field_shift()
    end

    function field_shifts()
    end
end
export info
