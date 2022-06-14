module operation

    export func, correct, significant
    using DataFrames, DataFramesMeta

    """
    func

    Inputs
    ------

    main, shuffle :: DataFrame

    transform :: Pair
        Which column to apply (main,shuffle) function to and what to call the
        new column

    shuffle_is_by :: set of df columns
        which dimensions should be in the shuffle
    func_dims :: set of df columns
        dimensions that should be in the final result
        (a union of shuffle_dims and these)

    """
    function func(main, shuffle, transform::Union{Symbol,Pair}; 
            func_dims=[], shuffle_is_by=[:shuffle,:unit], metric=:info,
            stat=:median, op=.-)
            
        stat = stat      isa Symbol ? eval(timeshift, stat) : stat
        tranform = transform isa Symbol ? (transform => transform) : transform
        inputfields, outputfields = transform

        @infiltrate
        shuffle = combine(groupby(shuffle, shuffle_is_by), 
                          inputfields => stat => inputfields)

        shuffle_is_by = [b for b in shuffle_is_by if b ∈ propertynames(main)]
        by = union(shuffle_is_by, func_dims)
        main, shuffle = groupby(main, by), groupby(shuffle, by)
        @assert size(main) == size(shuffle)

        result = DataFrame()
        for (m, s) in zip(main, shuffle)
            r = m[!, Not(statopfields)]
            r[!,outputfields] = op(m[!,inputfields], s[!,inputfields])
            append!(result, r)
        end

    end
    function significant(main, shuffle, transform::Union{Symbol, Pair}; kws...)
        kws = (;kws..., op = ((a,b) -> mean(a .> b , dims=1)))
        func(main, shuffle, transform; kws...)
    end
    function correction(main, shuffle, transform; stat=:mean, kws...)
        kws = (;kws..., op=.-, stat=stat)
        func(main, shuffle, transform; kws...)
    end

end
