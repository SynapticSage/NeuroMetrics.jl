module GoalFetchAnalysis

    __revise_mode__ = :eval

    # General
    include("raw.jl")   # Loading/saving/manipulating raw data
    include("utils.jl") # General utilites
    include("filt.jl")
    include("shuffle.jl") # General receptive field codes
    import .raw
    import .utils
    import .filt
    import .shuffle

    include("table.jl") # Manipulating tables
    import .table

    # Field related
    include("field.jl") # General receptive field codes
    import .field

    include("timeshift.jl")
    import .timeshift

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
