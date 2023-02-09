"""
# Plots for manifold related objects

#### See
maniplot
"""
module manifold

    using Plots
    using Infiltrator
    using RecipesBase
    using Statistics
    import Random
    using DIutils
    import DIutils.plotutils: getplotcolor

    # KEEP THIS HERE FOR INTRA-SCRIPT TESTING @userplot ManiPlot
    @recipe function maniplot(plt::ManiPlot; n=14_000, kalpha=0.5, 
            llim=nothing, origin=0, 
            by=nothing, byaxis=false,
            quantfilt=0.999,
            linealpha=nothing, pointalpha=nothing, cmap=:vik)

        em = plt.args[1]
        n = min(n, size(em,1))
        @assert em isa AbstractArray
        if quantfilt !== nothing
            qs = [ e .> quantile(e, 1-quantfilt) .&& e .< quantile(e, quantfilt)
                  for e in eachcol(E) ]
            qs = accumulate!(.&&, qs)
            E = E[qs, :]
            by = by[qs, :]
        end
        I = 1:n
        E = em[I,:];
        E = if byaxis 
            B = by[I]
            B = if llim === nothing
                B
            else
                B = DIutils.nannorm_extrema(B)
                B = llim*2(B) - origin
            end
            [E B] 
        else
            E
        end

        @series begin
            colormap --> cmap
            width := 10
            subplot := 2
            seriestype := :heatmap
            (fill(NaN,2,2))
        end
        if llim !== nothing
            @info llim
            ylims := (origin-llim, origin+llim)
            xlims := (origin-llim, origin+llim)
            zlims := (origin-llim, origin+llim)
        end
        @series begin 
            subplot --> 1
            if linealpha !== nothing
                alpha := linealpha
            end
            seriestype := :line
            color := :black
            label := nothing
            (eachcol(E)...,)
        end
        subplot --> 1
        xlabel --> "x"
        ylabel --> "y"
        zlabel --> "z"
        colorbar --> true
        alpha --> 0.1 * kalpha
        markersize --> 2
        seriestype := :scatter
        if pointalpha !== nothing
            alpha := pointalpha
        end
        color --> if by !== nothing
            @info cmap
            getplotcolor(by[I], cmap)
        else
            theme_palette(:auto)[1]
        end
        label := nothing
        (eachcol(E)...,)
    end

end
