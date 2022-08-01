"""
`rf`
"""
module receptivefield

    using Field
    using Field: ReceptiveField, Grid, Occupancy
    using Field.adaptive: GridAdaptive
    import Utils
    using Plots, LaTeXStrings, Measures
    using Statistics
    using ProgressMeter
    using NaNStatistics
    import TextWrap
    using Infiltrator
    



    @recipe function plot_adaptiverf(field::ReceptiveField, val::Symbol=:rate;
            ztransform::Bool=false, mfunc::Function=nanmean,
            sfunc::Function=nanstd, title_width=40)
        seriestype --> :heatmap
        title --> TextWrap.wrap(string(field.metrics), width=title_width)
        X = [field.grid.centers[1]...]
        x --> X
        if length(field.grid.centers) > 1
            Y = [field.grid.centers[2]...]
            y --> Y
        end
        Z = getproperty(field, val)
        if ztransform
            Z = (Z .- mfunc(Z))./sfunc(Z)
            colorbar_title --> String(val) * " Z-transform"
        else
            colorbar_title --> String(val)
        end

        (X, Y, Z')
    end

    @recipe function plot_adaptiveocc(field::T where T<:Occupancy,
            val::Symbol=:prob)
        seriestype --> :heatmap
        colorbar_title --> String(val)
        seriestype --> :heatmap
        X = [field.grid.centers[1]...]
        x --> X
        if length(field.grid.centers) > 1
            Y = [field.grid.centers[2]...]
            y --> Y
        end
        z := if val==:prob
            (X, Y, reshape(getproperty(field, val), size(field.grid))')
        else
            (X, Y, getproperty(grid, val)')
        end
    end

    @recipe function plot_adaptivegrid(grid::GridAdaptive, val::Symbol=:radii)
        colorbar_title --> String(val)
        seriestype --> :heatmap
        c --> :thermal
        X = [grid.centers[1]...]
        x --> X
        if length(grid.centers) > 1
            Y = [grid.centers[2]...]
            y --> Y
        end
        Z = getproperty(grid, val)
        if eltype(Z) == Vector
            Z = reshape([sqrt(sum(z.^2)) for z in Z], size(Z))
        end
        (X, Y, Z')
    end


end

