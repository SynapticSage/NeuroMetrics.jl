# fetch-maze-analysis-julia

## What's here?
64-tetrode fetch maze analysis. Anything I do with the julia language to analyze my high-denisty recordings. Separate from the matlab repo for the same dataset as of now.

This code base uses [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> goal-code

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all of the correct Julia libraries (used in plotting/munging)

2. Follow instructions for [[https://github.com/JuliaInterop/MATLAB.jl][MATLAB.jl]]

    you may need to set your Julia environment variable
    ```
    env["MATLAB_ROOT"] = "/usr/local/MATLAB/R2021b/";
    ```
    or wherever you store Matlab

