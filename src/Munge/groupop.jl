__precompile__(false)
module groupop

    export Distribute
    struct Distribute
        X
    end

    """
        function prod(f, d1::Distribute, d2::Distribute)
    Returns the product of the elements of d1 and d2, where
    f is a function that takes two arguments.
    # Arguments
    - `f`: Function that takes two arguments.
    - `d1::Distribute`: Distribute object.
    - `d2::Distribute`: Distribute object.
    # Returns
    - `Array`: Array of the products of the elements of d1 and d2.
    """
    function Base.prod(f, d1::Distribute, d2::Distribute)
        indices = collect(Iterators.product(eachindex(d1.X), eachindex(d2.X)))
        typ = typeof(f(d1.X[1], d2.X[1]))
        out = Array{typ}(undef, size(indices))
        Threads.@threads for (i, j) in indices
            out[i, j] = f(d1.X[i], d2.X[j])
        end
        return out
    end

    function _index(X)
        eachindex(X)
    end
    function _index(X::AbstractDict)
        keys(X)
    end
    function _index(X::GroupedDataFrame)
        keys(X)
    end


end
