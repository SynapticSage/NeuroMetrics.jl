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

    # KEEP THISHERE FOR INTRA-SCRIPT TESTING @userplot ManiPlot
    @userplot ManiPlot
    @recipe function maniplot(plt::ManiPlot; n=14_000, kalpha=0.5, 
            llim=nothing, origin=nothing, 
            by=nothing, byaxis=false,
            quantfilt=0.99,
            linealpha=nothing, pointalpha=nothing, cmap=:vik)

        em = plt.args[1]
        K = 1
        @assert em isa AbstractArray
        if quantfilt !== nothing
            qs = [ (e .> quantile(e, 1-quantfilt) .&& 
                    e .< quantile(e, quantfilt))
                  for e in eachcol(em) ]
            qs = hcat(qs...)
            @infiltrate
            qs = findall(vec(all(qs, dims=2)))
            em = em[qs, :]
            by = by !== nothing ? by[qs, :] : nothing
            #@info "size(em) = $(size(em))"
        end
        if origin !== nothing
            origin = length(origin) != size(em,2) ?
                fill(origin, size(em,2)) : origin
        else
            origin = median(em, dims=1)
            @info "median, $origin"
        end
        if llim !== nothing
            llim = length(llim) != size(em,2) ?
                fill(llim, size(em,2)) : llim
        end

        n = min(n, size(em,1))
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
        if llim !== nothing
            #@info llim
            ylims := (origin[1]-llim[1], origin[1]+llim[1])
            xlims := (origin[2]-llim[2], origin[2]+llim[2])
            zlims := (origin[3]-llim[3], origin[3]+llim[3])
        end
        @series begin 
            subplot --> 1
            if linealpha !== nothing
                alpha := linealpha
            else
                alpha := kalpha * K
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
        alpha --> kalpha * K
        markersize --> 2
        seriestype := :scatter
        if pointalpha !== nothing
            alpha := pointalpha
        end
        color --> if by !== nothing
            #@info cmap
            getplotcolor(by[I], cmap)
        else
            theme_palette(:auto)[1]
        end
        label := nothing
        (eachcol(E)...,)
    end

end
