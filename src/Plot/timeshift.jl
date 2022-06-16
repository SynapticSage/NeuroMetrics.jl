module timeshift

    using StatsPlots: @df
    using Plots
    using Infiltrator
    import Shuffle
    using DataFrames, DataFramesMeta

    export plot_shifts

    dosave = true

    function _xlim(df::DataFrame)
        if :shift ∈ propertynames(df)
            xlim = (nanminimum(df.shift),   nanmaximum(df.shift))
        elseif :bestTau ∈ propertynames(df)
            xlim = (nanminimum(df.bestTau), nanmaximum(df.bestTau))
        else
            throw(KeyError("Missing valid column for auto-xlim in dataframe"))
        end
        return xlim
    end
    xlabel = "shift\n(seconds)"

    #function _plot(df::DataFrame, pos::Symbol...; kws...)
    #    xlim = get(kws, :xlim, _xlim(df_imax))
    #    seriestype = get(kws, :seriestype, :density)
    #    @df df plot(pos...; 
    #end

    function plot_shuffle(df::DataFrame, method; shuffle=true, takefirst=10, display=:rows, kws...)
        if method isa Symbol
            method = eval(method)
        end
        df = collect(groupby(df, :shuffle))
        if shuffle
            Shuffle.shuffle!(df)
        end
        plots = []
        for d in Iterators.take(df, takefirst)
            push!(plots,method(DataFrame(d); kws...))
        end

        @assert length(plots) == takefirst

        if display==:overlay
            P = Plots.plot(pop!(plots))
            for p in plots
                Plots.plot!(p)
            end
            P
        elseif display==:rows
            Plots.plot(plots...; grid=grid(takefirst,1))
        elseif display==:cols
            Plots.plot(plots...; grid=grid(takefirst,1))
        else
            throw(ArgumentError("Unrecognized display=$display"))
        end
    end

    function plot_bestTauEst(df_imax::DataFrame; desc::String="", kws...)
        xlim = get(kws, :xlim, _xlim(df_imax))
        seriestype = get(kws, :seriestype, :density)
        group=get(kws, :group, df_imax.area)
        @df df_imax plot(:bestTau; 
                            seriestype=seriestype,
                            title="$desc BestTau(Information)", 
                            xlim, xlabel, ylabel="Density", kws...)
    end

    function plot_infoMean(df::DataFrame; desc::String="", kws...)
        xlim = get(kws, :xlim, _xlim(df))
        seriestype = get(kws, :seriestype, :density)
        if :info_mean ∉ propertynames(df)
            @info "mean-ing"
            df = info_mean(df)
        end
        group=get(kws, :group, df.area)
        @df df plot(:shift, :info_mean; group, xlim, 
                    title="$desc Median(Information)", xlabel, 
                    ylabel="Shannon MI Bits", kws...)
    end

    function plot_distributionPerCell(df::DataFrame; measure=:info, 
            statistic=mean, sortorder=nothing,
            clim=:cell, desc::String="")

        @infiltrate
        df = sort(df, [:shift,:area,:unit])
        if :shuffle in propertynames(df) && length(unique(df.shuffle))>1
            start = 4
            DF = combine(groupby(df, [:unit,:shift]), 
                         Not([:unit,:shift,measure]) .=> first .=> Not([:unit,:shift,measure]), 
                         measure => statistic => measure)
        else
            DF = df
            start = 3
        end
        df_u = sort(unstack(DF, :shift, measure), [:area, :unit])
        units = df_u.unit
        areas = df_u.area
        area_divide = findfirst([diff(areas .== "PFC"); 0].==1)
        valcols = [x for x in names(df_u) if Symbol(x) ∉ [:area,:unit,:taus, :shuffle]]
        df_u.taus = vec([x[2] for x in argmax(Matrix(df_u[!,valcols]),dims=2)])
        if sortorder != nothing
            df_u = sort(df_u, [:area, :taus, :unit])
        else
            df_u = sort(df_u, [:area, :taus, :unit])
        end
        shifts = parse.(Float32, valcols)

        function get_area(area) 
            # cells x shifts
            M = Matrix(df_u[df_u.area.==area, valcols])
            M = replace(M, missing=>NaN)
            M[findall(vec([x[1] for x in any(M.!=0,dims=2)])),:]
        end

        norm₁(x, m, M) = begin
            #@infiltrate nanmaximum(x) > 0
            (max.(min.(x,M),m) .- m)./(M-m)
        end

        ca1, pfc = get_area("CA1"), get_area("PFC")

        @assert length(shifts) == size(ca1,2)
        if clim == :area
            ca1_clim = nanquantile.(vec(ca1), [0.01, 0.95])
            pfc_clim = nanquantile.(vec(pfc), [0.01, 0.95])
        elseif clim == :cell
            ca1_clim = hcat(nanquantile(ca1, 0.01, dims=2),
                            nanquantile(ca1, 0.95, dims=2))
            pfc_clim = hcat(nanquantile(pfc, 0.01, dims=2),
                            nanquantile(pfc, 0.95, dims=2))
            for (i,cc) in zip(1:size(ca1,1), eachrow(ca1_clim))
                ca1[i,:] = norm₁(ca1[i,:], cc[1], cc[2])
            end
            for (i, pp) in zip(1:size(pfc,1), eachrow(pfc_clim))
                pfc[i,:] = norm₁(pfc[i,:], pp[1], pp[2])
            end
            ca1_clim = (0,1)
            pfc_clim = (0,1)
        else
        end

        p = Plots.plot(
                   heatmap(shifts, 1:size(get_area("CA1"),1), ca1, 
                    clims=ca1_clim, colorbar_title="Spatial-firing mutual information", colorbar_titlefontrotation=0),
                   heatmap(shifts, 1:size(get_area("PFC"),1), pfc, 
                    clims=pfc_clim, colorbar_title="Spatial-firing mutual information", colorbar_titlefontrotation=0),
                   title="$desc\nMI", xlabel=xlabel, ylabel="cell"
        )
        vline!(p[1], [0], c=:white, linestyle=:dash, label="Zero lag")
        vline!(p[2], [0], c=:white, linestyle=:dash, legendposition=:none)
    end
    
    function saveplotshift(args...; formats=["png","svg","pdf"])
        args = ["fields", "shifts", args...]
        if dosave
            for format in formats
                savefig(plotsdir(args[1:end-1]..., args[end] * ".$format"))
            end
        end
    end

    # --------
    # PLOTTING
    # --------
    function plot_shifts(place::AbstractDict; desc::String="", shift_scale=:minutes, clim=:cell)

        # List of desirables
        # ------------------
        # - shuffle plots
        #   - all of the below plots with shuffle
        #       - correction
        #       - probability
        # - unify time annotations across plots
        # - normalized version of plots :: where information is [0:lowest for neuron,
        #                                                       1:highest for neuron]
        # - quantile based clims
        #   currently clipping may be changing my gut feeling about the data
        # - save the imax dataframe

        descSave = replace(desc, ":"=>"", " "=>"-")

        # Setup our datatypes
        df = info_to_dataframe(place; shift_scale)
        df_imax  = imax(df)

        # Some settings across everything

        # Density of bestTau point estimates per cell
        density_bestTauEst(df_imax; desc)
        saveplotshift("$(descSave)_density_x=shift,y=densBestTau_by=area")

        # HISTOGRAM ALL CELLS
        df_m = info_mean(df)
        saveplotshift("$(descSave)_histogram_x=shift,y=medianinfo_by=area")

        # PER CELL
        saveplotshift("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.pdf")
    end

end
