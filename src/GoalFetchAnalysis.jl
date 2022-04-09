module GoalFetchAnalysis

    include("raw.jl")   # Loading/saving/manipulating raw data
    include("table.jl") # Manipulating tables
    include("raster.jl")# Plotting rasters
    include("field.jl") # General receptive field codes
    include("utils.jl") # General utilites
    include("field/operation.jl") # Operations/munging on receptive fields

    export raw
    export table 
    export raster 
    export field

end
