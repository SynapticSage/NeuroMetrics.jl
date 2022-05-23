
module group
    using DataFrames
    using DataStructures
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
    function multitable_groupby(groups::Union{Vector, Symbol, String}, args::AbstractDataFrame...)
        args = vcat(args...; source="source")
        [groupby(g, :source) for g in groupby(args, groups)]
    end
end
export group

