module timeshift

    using StatsPlots: @df
    using Plots
    using Infiltrator
    import Shuffle
    import Table
    import Utils
    using Statistics, NaNStatistics
    using DataFrames, DataFramesMeta
    using StatsBase

    export plot_shifts
    export plot_shift_versus_info

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

    function heatmap_unitshift(Isc::DataFrame; sort_rows=:bestshift,
            value=:value,
            removenonsig=false, centerline=true, setnanzero=true,
            dropmissingrows=false, kws...)
		 if !(typeof(sort_rows) <:Vector)
					sort_rows = [sort_rows]
		end
        extra = removenonsig ? [:sig_minimum] : []
		not_list = union([:unit, extra...], sort_rows)
		shifts = unique(sort(Isc.shift))
		dropmissing!(Isc)
		ustack = unstack(sort(Isc,[:unit,:shift, extra...]), not_list, :shift, value)
		ustack = sort(ustack, sort_rows)
        if removenonsig
            ustack = @subset(ustack, :sig_minimum .< 0.05)
        end
		units  = unique(sort(ustack.unit))
		ustack = ustack[:,Not(not_list)]
		ustack = replace(Matrix(ustack), missing=>NaN)	
        good = (!).(isnan.(ustack))
        if dropmissingrows && !(setnanzero)
            inds= findall(Utils.squeeze(any(good; dims=2)))
            ustack = ustack[inds, :]
            units = units[inds]
        elseif setnanzero
            bad = findall(Utils.squeeze((!).(good)))
            ustack[bad,:] .= 0
            #units = units[inds]
        end
        hm = Float32.(Matrix(ustack))
        p = heatmap(shifts, 1:length(units), hm; kws...)
        if centerline
            vline!([0], c=:white, linestyle=:dash, linewidth=2, label="")
        end
        p
    end

    """
    Calculates the empirical cdf of the cell relative to its shuffle
    """
    function ecdf_units(I, S; groups=:unit)

        G = Table.group.mtg_via_commonmap(groups, I, S)
        for (ig, sg) in G
            # Plot
            C = ecdf(sg.value)
            p =plot(C, title="unit=$(Int64(ig.unit[1])) area=$(ig.area[1])",
                                xlabel="bits per spike", ylabel="fraction", ylim=[0,1])
            plot!([ig.value..., ig.value...], [ylims()...])
            annotation = string(round(Float64(ig.value[1]),digits=2))
            @infiltrate annotation != "NaN"
            annotate!(ig.value, mean([ylims()...]), text(annotation))
            if ig.sig[1] < 0.05
                push!(P_sig, p)
            else
                push!(P_non, p)
            end
        end

        P_sig, P_non

    end

    function scatter_unitfield_by_unitfield(Is, x, y; kws...)
        groups = combine(groupby(Is, :unit), first)
        s = sort(groups[:,[:unit, y, x]], y)
        scatter(s[:,x], s[:,y]; xlabel=x, ylabel=y, kws...)
    end

    function plot_shift_versus_info(Is, sig_frac, bonf_cells_are_sig; 
                                     groups=:shift, kws...)
        subset_of_interest = subset(bonf_cells_are_sig, :sig => x -> x .== sig_frac)
        sig_cells = Utils.ismember(Is.unit, subset_of_interest.unit)
        #println("Insig cell count  = $(round(sum(insig_cells)/length(unique(Is.shift))))")
        sig_cells = Is[sig_cells, :]
        # Candidate for a Util
        sig_shape = combine(groupby(sig_cells, groups), 
                :value => mean,
                :value => median,
                :value => (x->nanquantile(x, 0.05))  => :lower,
                :value => (x->nanquantile(x, 0.95))  => :upper)
        sort!(sig_shape, :shift)
        p = @df sig_shape plot(:shift * 60, :value_mean; 
                               fillrange=(:lower, :upper), label="mean", kws...)
        @df sig_shape plot!(:shift*60, :lower, style=:dash, c=:black, label="lower")
        @df sig_shape plot!(:shift*60, :upper, style=:dash, c=:black, label="upper")
        hline!([0], c=:black, label="")
        p
    end

end
