module group

    using DataFrames
    using DataStructures
    using Infiltrator
    coords(groups) = collect(zip(sort(collect(groups.keymap),by=x->x[2])...))[1]
    pairs(groups) = (collect(zip(sort(collect(groups.keymap))))[i][1] 
                     for i in 1:length(G.keymap))
    function named_coords(groups)
    end

    """
    finds a common index system to move through multiple dataframes
    """
    function find_common_mapping(args::GroupedDataFrame...)
        # Obtain a hash of the keys of each keymap
        hashmaps = []
        for i ∈ 1:length(args)
            hashmap = OrderedDict()
            for (key, value) in args[i].keymap
                hashmap[hash(key)] = value
            end
            push!(hashmaps, hashmap)
        end
        # Intersect the hashmaps
        superset = [keys(hashmap) for hashmap in hashmaps]
        superset = Tuple(intersect(superset...))
        # Throw out keys that fail to match the superset
        for i ∈ 1:length(args)
            hashmap = OrderedDict(key=>hashmaps[i][key] for key in superset)
            hashmaps[i] = hashmap
        end
        # Now that all of the orderings match, I can read out the matrix of
        # group numbers per matching key
        result = [hashmaps[i][superset[j]] for j ∈ 1:length(superset), i ∈ 1:length(args)]
    end

    """
    Combine two tables and groupby some props, while distinguishing the tables with a new key
    """
    function multitable_groupby(groups::Union{Vector, Symbol, String}, args::AbstractDataFrame...; dropmissingrows::Bool=false)
        args = vcat(args...; source="source", cols=:union)
        G = [groupby(g, :source) for g in groupby(args, groups)]
        source_count, group_count = maximum([length(G[i]) for i in 1:length(G)]),
                                    length(G)
        GG = Array{Union{Missing, SubDataFrame}}(undef, group_count, source_count)
        for i in 1:group_count
            for j in 1:length(G[i])
                GG[i,j] = G[i][j]
            end
            if length(G[i]) < source_count
                for j in (length(G[i])+1):source_count
                    GG[i,j] = missing
                end
            end
        end
        if dropmissingrows
            GG = GG[sum(ismissing.(G),dims=2).==0,:]
        end
        GG
    end

end

