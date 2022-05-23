

    # --------
    # PLOTTING
    # --------
    function plot_shifts(place; desc="", shift_scale=:minutes, clim=:cell)

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
        function saveplotshift(args...)
            savefig(plotsdir(args[1:end-1]..., args[end] * ".png"))
            savefig(plotsdir(args[1:end-1]..., args[end] * ".svg"))
            savefig(plotsdir(args[1:end-1]..., args[end] * ".pdf"))
        end

        # DENSITY ALL CELLS
        df, df_imax = info_dataframe_and_cell_dataframe(place; shift_scale)
        xlabel = "shift\n(seconds)"
        xlim = (nanminimum(df.shift), nanmaximum(df.shift))

        @df df_imax density(:bestTau, group=:area, 
                     title="$desc BestTau(Information)", xlim=xlim,
                     xlabel=xlabel, ylabel="Density")
        saveplotshift("fields", "shifts",
             "$(descSave)_density_x=shift,y=densBestTau_by=area")
        @df df_imax histogram(:bestTau; group=:area,  xlim,
                     title="$desc BestTau(Information)", 
                     xlabel=xlabel, ylabel="Density")
        saveplotshift("fields", "shifts",
             "$(descSave)_histogram_x=shift,y=densBestTau_by=area")

        # HISTOGRAM ALL CELLS
        df_m = combine(groupby(df, [:area, :shift]), :info=>mean)
        @df df_m bar(:shift, :info_mean; group=:area, xlim,
                     title="$desc Median(Information)", 
                     xlabel=xlabel, ylabel="Shannon MI Bits")
        saveplotshift("fields", "shifts",
             "$(descSave)_histogram_x=shift,y=medianinfo_by=area.svg")

        # PER CELL
        df = sort(df, [:shift,:area,:unit])
        df_u = sort(unstack(df, :shift, :info), [:area, :unit])
        shifts = parse.(Float32,names(df_u)[3:end])
        units = df_u.unit
        areas = df_u.area
        area_divide = findfirst([diff(areas .== "PFC"); 0].==1)
        exclude = Not([x for x in [:area,:unit,:taus]
                       if String(x) in names(df_u)])
        df_u.taus = vec([x[2] for x in argmax(Matrix(df_u[!,exclude]),dims=2)])
        df_u = sort(df_u, [:area, :taus])
        function get_area(area) 
            # cells x shifts
            M = Matrix(df_u[df_u.area.==area, Not([:area,:unit, :taus])])
            M = replace(M, missing=>NaN)
            M[findall(vec([x[1] for x in any(M.!=0,dims=2)])),:]
        end
        norm₁(x, m, M) = begin
            #@infiltrate nanmaximum(x) > 0
            (max.(min.(x,M),m) .- m)./(M-m)
        end
        ca1, pfc = get_area("CA1"), get_area("PFC")
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
        saveplotshift("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.pdf")
    end
