module GoalFetchAnalysis

    include("raw.jl")
    include("table.jl")
    include("raster.jl")
    include("field.jl")
    include("utils.jl")

    export raw
    export table 
    export raster 
    export field

end
