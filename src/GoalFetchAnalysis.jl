module GoalFetchAnalysis

    import Distributed: @everywhere
    @everywhere using DrWatson
    #push!(LOAD_PATH, @__DIR__)
    __revise_mode__ = :eval
    @everywhere push!(LOAD_PATH, srcdir())

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

    import Utils
    import Table
    import Load, Filt, Shuf
    import Field
    #import Decode
    import Munge
    import Timeshift

    export Utils
    export Table 
    export Filt
    export Load
    export Field
    export Timeshift
    export Decode
    export Munge
    export Plot

    filtreg = Utils.filtreg

    #Utils.plot.set_theme_timebased(23)
    import Plots
    Plots.theme(:bright)

    function import_timeshift()
        @eval Main using Timeshift
        @eval Main using Timeshift.checkpoint
    end

    import Labels as labels

    # EVERY time I import this, I check my tests to ensure
    # I've broken nothing during development
    # Run quick tests at startup to ensure library integrity
    if isdir(projectdir("test")) && isfile(projectdir("test","runtests.jl")) && 
        "USE_PLUTO" ∉ keys(ENV)
        include(projectdir("test","runtests.jl"))
    end
    
end
