module GoalFetchAnalysis

    using Revise, DrWatson
    __revise_mode__ = :eval
    import Plots
    Plots.theme(:bright)

    import DIutils 
    import DIutils.Table as Table
    export DIutils, Table

    import DI
    import DI: Filt, Labels
    export DI, Filt

    include("Shuf.jl")
    include("Field.jl")
    include("Timeshift.jl")
    include("Munge.jl")
    include("Plot.jl")
    # incldue("Decode.jl")

    export Filt, Labels
    export DI
    export Field
    export Timeshift
    export Munge
    export Plot
    # export Decode

    filtreg = DIutils.filtreg

    #import Labels as labels
    labels = DI.Labels

    # ------
    # REVISE.jl
    # ------
    # EVERY time I import this, I check my tests to ensure
    # I've broken nothing during development
    # Run quick tests at startup to ensure library integrity
    function check_dir_is_watched(dir)
        dir = srcdir(dir)
        watched = keys(Revise.watched_files)
        if endswith(dir, '/')
            dir = dir[1:end-1]
        end
        dir_in_watched = dir ∈ watched
        !(dir_in_watched) ? @warn("dir=$dir not in watched") : nothing
        dir_in_watched
    end
    export pushover
    pushover = DIutils.pushover
    pushover("GoalFetchAnalysis.jl ready")

    include("Precompile.jl")


end

    # import Distributed: @everywhere
    #push!(LOAD_PATH, @__DIR__)
    #@everywhere push!(LOAD_PATH, srcdir())

    ## General
    #include(srcdir("DI.jl"))
    #using .DI   # Loading/saving/manipulating raw data
    #include(srcdir("DIutils.jl"))
    #using .DIutils # General utilites
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

    
    #if isdir(projectdir("test")) &&isfile(projectdir("test","runtests.jl")) && 
    #   "USE_PLUTO" ∉ keys(ENV)
    #   include(projectdir("test","runtests.jl"))
    #end
    # DIutils.plot.set_theme_timebased(23)

    # function import_timeshift()
    #     @eval Main using Timeshift
    #     @eval Main using Timeshift.checkpoint
    # end

