module GoalFetchAnalysis

    push!(LOAD_PATH, @__DIR__)
    module_dir(d) = isdir(d) && isuppercase(d[1])
    push!(LOAD_PATH, joinpath(@__DIR__, "Load"))
    push!(LOAD_PATH, joinpath(@__DIR__, "Utils"))
    push!(LOAD_PATH, joinpath(@__DIR__, "Table"))
    push!(LOAD_PATH, joinpath(@__DIR__, "Field"))
    push!(LOAD_PATH, joinpath(@__DIR__, "Filt"))
    push!(LOAD_PATH, joinpath(@__DIR__, "Shuffle"))
    @info LOAD_PATH

    __revise_mode__ = :eval
    using OhMyREPL

    # General
    using Load   # Loading/saving/manipulating raw data
    using Utils # General utilites
    import Filt
    import Shuffle # General receptive field codes

    import Table # Manipulating tables

    # Field related
    import Field # General receptive field codes
    import Timeshift

    # Goal vector measures
    import statistic

    # Deode Related
    import Decode # General receptive field codes

    export Load
    export Utils
    export Table 
    export Field, Timeshift
    export Decode

end
