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
            sfunc::Function=nanstd, title_width=40, transpose::Bool=true)
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

        transpose ? (Y, X, Z) : _transpose(Y,X,Z)
    end

    @recipe function plot_adaptiveocc(field::T where T<:Occupancy,
            val::Symbol=:prob, transpose::Bool=true)
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
            (Y, X, getproperty(grid, val))
        end
        transpose ? (X,Y,Z') : _transpose(X,Y,Z')
    end

    @recipe function plot_adaptivegrid(grid::GridAdaptive, val::Symbol=:radii,
        transpose::Bool=true)
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
        transpose ? (X,Y,Z') : _transpose(X,Y,Z')
    end

    function _transpose(X, Y, Z)
        (Y,X,Z')
    end


end

