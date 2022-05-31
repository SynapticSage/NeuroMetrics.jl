module GoalFetchAnalysis

    # General
    include("raw.jl")   # Loading/saving/manipulating raw data
    include("table.jl") # Manipulating tables
function cellpath(animal::String, day::Int, tag::String=""; type="csv", kws...)
    if tag != "" && tag != "*"
        tag = "_$tag"
    end
    path = DrWatson.datadir("exp_raw", "visualize_raw_neural",
                               "$(animal)_$(day)_cell$tag.$type")
    path = cellpath(pos...; kws...)
    if occursin("*", path)
        base, dir = basename(path), dirname(path)
        @debug "base=$base, dir=$dir"
        paths = glob(base, dir)
    else
        paths = [path]
    end
end

    include("utils.jl") # General utilites
    include("filt.jl")
    include("shuffle.jl") # General receptive field codes
    import .raw
    import .table
    import .field
    import .utils
    import .filt
    import .shuffle

    # Field related
    include("field.jl") # General receptive field codes
    import .field
    timeshift = field.timeshift

    # Goal vector measures
    include("statistic.jl")
    import .statistic

    # Deode Related
    include("decode.jl") # General receptive field codes
    import .decode

    using OhMyREPL

    export raw
    export table 
    export field, timeshift
    export decode
    export utils

end
