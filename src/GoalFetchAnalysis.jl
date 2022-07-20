module GoalFetchAnalysis

    using DrWatson
    #push!(LOAD_PATH, @__DIR__)
    __revise_mode__ = :eval
    push!(LOAD_PATH, srcdir())

    ## General
    #include(srcdir("Load.jl"))
    #using .Load   # Loading/saving/manipulating raw data
    #include(srcdir("Utils.jl"))
    #using .Utils # General utilites
    #include(srcdir("Filt.jl"))
    #import .Filt
    #include(srcdir("Shuffle.jl"))
    #import .Shuffle # General receptive field codes
    #include(srcdir("Table.jl"))
    #import .Table # Manipulating tables
    ## Field related
    #include(srcdir("Field.jl"))
    #import .Field # General receptive field codes
    #include(srcdir("Timeshift.jl"))
    #import .Timeshift
    ## Goal vector measures
    #include(srcdir("Statistic.jl"))
    #import .Statistic
    ## Deode Related
    #include(srcdir("Decode.jl"))
    #import .Decode # General receptive field codes
    #include(srcdir("Plot.jl"))
    #import .Plot

    import Load, Filt, Shuf, Table, Field
    import Timeshift
    import Decode, Munge, Utils

    export Load
    export Utils
    export Table 
    export Field
    export Timeshift
    export Decode
    export Filt
    export Plot
    export Munge

    Utils.plot.set_theme_timebased(23)

end
