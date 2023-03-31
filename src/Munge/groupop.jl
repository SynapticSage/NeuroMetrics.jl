__precompile__(false)
module groupop

    using DataFrames
    using Infiltrator

    export AllToAll
    """
        struct AllToAll
    A struct that allows all elements of a collection to be
    operated on by all elemnets of a another collection
    """
    struct AllToAll
        X
    end
    Base.getindex(d::AllToAll, i) = d.X[i]
    Base.setindex!(d::AllToAll, v, i) = d.X[i] = v
    Base.eachindex(d::AllToAll) = eachindex(d.X)
    Base.axes(d::AllToAll) = axes(d.X)

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
    function Base.map(f, d1::AllToAll, d2::AllToAll)
        indices = collect(Iterators.product(eachindex(d1.X), eachindex(d2.X)))
        typ = typeof(f(d1.X[1], d2.X[1]))
        out = Array{typ}(undef, size(indices))
        Threads.@threads for (i, j) in indices
            out[i, j] = f(d1.X[i], d2.X[j])
        end
        return out
    end

    Base.:*(d1::AllToAll, d2::AllToAll) = map(.*, d1, d2)
    Base.:+(d1::AllToAll, d2::AllToAll) = map(.+, d1, d2)
    Base.:-(d1::AllToAll, d2::AllToAll) = map(.-, d1, d2)
    Base.:/(d1::AllToAll, d2::AllToAll) = map(./, d1, d2)

    struct OneToOne
        X
    end
    Base.getindex(d::OneToOne, i) = d.X[i]
    Base.setindex!(d::OneToOne, v, i) = d.X[i] = v
    Base.eachindex(d::OneToOne) = eachindex(d.X)
    Base.axes(d::OneToOne) = axes(d.X)
    Base.size(d::OneToOne) = size(d.X)
    Base.length(d::OneToOne) = length(d.X)

    function Base.map(f, d1::OneToOne, d2::OneToOne)
        out = Array{typ}(undef, size(indices))
        Threads.@threads for (i, j) in eachindex(d1.X)
            out[i, j] = f(d1.X[i], d2.X[j])
        end
        return out
    end

    Base.:*(d1::OneToOne, d2::OneToOne) = map(*, d1, d2)
    Base.:+(d1::OneToOne, d2::OneToOne) = map(+, d1, d2)
    Base.:-(d1::OneToOne, d2::OneToOne) = map(-, d1, d2)
    Base.:/(d1::OneToOne, d2::OneToOne) = map(/, d1, d2)


    # Conversions
    AllToAll(X::OneToOne) = AllToAll(X.X)
    OneToOne(X::AllToAll) = OneToOne(X.X)

end
