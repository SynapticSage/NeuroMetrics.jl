module Precompile

    export precompile_goalmaze
    export precompile_GFA_dependencies
    # Uses https://github.com/JuliaLang/PackageCompiler.jl to speed up use of my library
    # Spends less time in precompile when a function/module is first called

    _GFA_dependencies=["DataFrames","Plots", "Makie", "CairoMakie", "OhMyREPL",
                       "CSV", "HDF5", "NetCDF", "Arrow", "Revise", "Reexport",
                       "DataFramesMeta", "Statistics", "NaNStatistics",
                       "ProgressMeter", "Glob", "Printf", 
                       "StatsPlots", "StatsBase", "Distributions", 
                       "Shuffle", "VideoIO",
                       "DataStructures", "Gadfly", "TableView", "LazyGrids",
                       "GeometricalPredicates", "RecipesBase", "Reexport",
                       "AxisArrays", "DimensionalData", "HypothesisTests",
                       "LoopVectorization", "Infiltrator",
                       "Markdown", "InteractiveUtils", "Polyester",
                       "Random", "Colors", "ColorSchemes", "ImageFiltering",
                       "ProtoStructs", "KernelDensity",
                       "UMAP", "Entropies", "CausalityTools"
                       "Term", # -- 27 jan 2023
                       "Images", "ImageSegmentation", "ImageTransformations", "ImageFiltering"
                      ]   
    _GFA_data_dependencies = [_GFA_dependencies..., "Load", "SampleData"]

    using PackageCompiler
    function precompile_GFA_dependencies(;incremental::Bool=true)
        PackageCompiler.create_sysimage(_GFA_dependencies;
                                        sysimage_path="GFA-dependencies-sysimage.so",
                                        incremental)
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
