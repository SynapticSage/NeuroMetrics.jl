# Uses https://github.com/JuliaLang/PackageCompiler.jl to speed up use of my library
# Spends less time in precompile when a function/module is first called

_GFA_dependencies=["DataFrames","Plots", "Makie", "GLMakie", "CairoMakie", "OhMyREPL",
                     "CSV", "HDF5", "NetCDF", "Arrow", "Revise",
                     "DataFramesMeta", "Statistics", "NaNStatistics",
                     "ProgressMeter", "Glob", "Printf", 
                     "StatsPlots", "StatsBase", "Distributions", "Shuffle",
                     "VideoIO", "DataStructures", "Gadfly", "Blink",
                     "ElectronDisplay", "TableView", "LazyGrids", "Random",
                     "Colors", "ColorSchemes", "ImageFiltering",
                     "KernelDensity"]   

using PackageCompiler
function precompile_GFA_dependencies()
    PackageCompiler.create_sysimage(_GFA_dependencies; sysimage_path="GFA-dependencies-sysimage.so")
end
function precomile_goalmaze()
    PackageCompiler.create_sysimage(["GoalFetchAnalysis"]; sysimage_path="GoalFetchAnalysis-sysimage.so")
end
