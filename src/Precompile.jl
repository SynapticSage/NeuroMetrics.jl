__precompile__(false)
module Precompile

    export precompile_goalmaze
    export precompile_GFA_dependencies
    import DIutils
    # Uses https://github.com/JuliaLang/PackageCompiler.jl to speed up use of my library
    # Spends less time in precompile when a function/module is first called

    _GFA_dependencies=["DataFrames", 
                        "FileIO",
                        #"UnicodePlots", # added unicodeplots for plots
                        "Plots", 
                       "Makie", "OhMyREPL", 
                       "CSV", 
                        # "HDF5", # removed 2023-05-01, check back later
                        "NetCDF", "Arrow", "Revise", 
                       #"Reexport","CairoMakie", 
                       "DataFramesMeta", "Statistics", "NaNStatistics",
                       "ProgressMeter", "Glob", "Printf", "StatsPlots",
                       "StatsBase", "Distributions", "Shuffle", "VideoIO",
                       "DataStructures", "Gadfly", "TableView", "LazyGrids",
                       "GeometricalPredicates", "RecipesBase",
                       "DimensionalData", "HypothesisTests",
                       "LoopVectorization", "Infiltrator", "Markdown",
                       "InteractiveUtils", "Polyester", "Random", "Colors",
                       "ColorSchemes", "ProtoStructs", "KernelDensity", # "UMAP",
                       "Entropies", "CausalityTools", "Term", "Images",
                       "ImageSegmentation", "ImageTransformations",
                       "ImageFiltering", "ArgParse", 
                       #"DIutils"
                        "MLJ" # 03/2023
]   
    _GFA_data_dependencies = [_GFA_dependencies..., "DI", "SampleData"]

    using PackageCompiler
    function precompile_GFA_dependencies(;incremental::Bool=true, append="")
        println("Precompiling GFA dependencies")
        println(_GFA_dependencies)
        PackageCompiler.create_sysimage(_GFA_dependencies;
                                        sysimage_path="GFA-dependencies-sysimage$append.so",
                                        incremental)
        DIutils.pushover("finished precompiling GFA dependencies at GFA-dependencies-sysimage$append.so")
    end
    function precompile_GFA_data_dependencies()
        PackageCompiler.create_sysimage(_GFA_data_dependencies;
                                        sysimage_path="GFA-dependencies-sysimage.so")
    end
    function precompile_goalmaze()
        PackageCompiler.create_sysimage(["GoalFetchAnalysis"];
                                        sysimage_path="GoalFetchAnalysis-sysimage.so")
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
