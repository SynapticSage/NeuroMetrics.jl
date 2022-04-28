module GoalFetchAnalysis

    # General
    include("raw.jl")   # Loading/saving/manipulating raw data
    include("table.jl") # Manipulating tables
    include("utils.jl") # General utilites
    incldue("filt.jl")
    include("shuffle.jl") # General receptive field codes
    import .raw
    import .table
    import .field
    import .utils
    import .filt

    # Field related
    include("field.jl") # General receptive field codes
    import .field

    # Goal vector measures
    include("statistic.jl")
    import .statistic


    export raw
    export table 
    export raster 
    export field

end
