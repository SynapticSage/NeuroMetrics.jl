module Precompile

    export precompile_goalmaze
    export precompile_GFA_dependencies
    # Uses https://github.com/JuliaLang/PackageCompiler.jl to speed up use of my library
    # Spends less time in precompile when a function/module is first called

    _GFA_dependencies=["DataFrames","Plots", "Makie", "CairoMakie", "OhMyREPL",
                       "CSV", "HDF5", "NetCDF", "Arrow", "Revise", "Reexport",
                       "DataFramesMeta", "Statistics", "NaNStatistics",
                       "ProgressMeter", "Glob", "Printf", "StatsPlots",
                       "StatsBase", "Distributions", "Shuffle", "VideoIO",
                       "DataStructures", "Gadfly", "TableView", "LazyGrids",
                       "Random", "Colors", "ColorSchemes", "ImageFiltering",
                       "ProtoStructs", "KernelDensity"]   

    using PackageCompiler
    function precompile_GFA_dependencies()
        PackageCompiler.create_sysimage(_GFA_dependencies; sysimage_path="GFA-dependencies-sysimage.so")
    end
    function precompile_goalmaze()
        PackageCompiler.create_sysimage(["GoalFetchAnalysis"]; sysimage_path="GoalFetchAnalysis-sysimage.so")
    end

    # -------------
    # Known issues
    # ------------
    #
    # precompile seems to defeat the normal use of electrondisplay, where once imported it takes over the display of plots and doc files.
    #
    # Things to try
    # - Remove electrondisplay from precompile
    # - Remove plots from precompile
end
