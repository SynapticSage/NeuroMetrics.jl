"""
Correlation(τ) and Correlation(rateᵢⱼ)

TODOs
----
- Higher resolution for XY analyses

"""

using DrWatson
using Timeshift
using Timeshift.types
using Timeshift.shiftmetrics
using Field.metrics
using Plot
using Plot.receptivefield
using Utils.namedtup
using DimensionalData
using ProgressMeter
import Plot
using Munge.spiking
using Filt
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics
using Plots
using GLM
using StatsBase, StatsPlots
using ElectronDisplay
using ProgressMeter
Plot.setfolder("timeshift","functional_connectivity")

@time spikes, beh, cells = Load.load("RY16", 36, data_source=["spikes","behavior", "cells"])
@time I = load_mains();
@time F = load_fields();
filt = get_filters()

metricfilters = Dict(
    :CA1=> x->
    x[:area] .=="CA1" .&& x[:meanrate] .> 0.01 && x[:meanrate] .< 7 .&& x[:maxrate] < 35 .&& x[:coherence] > 0.5,
    :PFC => x->
    x[:area] .=="PFC" .&& x[:meanrate] .> 0.005 .&& x[:maxrate] < 100,
   )
metricfilters[:CA1PFC] = x-> metricfilters[:CA1](x) || metricfilters[:PFC](x)

@showprogress for (datacut, region) in Iterators.product(
    [:all, :cue, :memory, :task, :nontask, :cue_correct, :cue_error, :mem_correct, :mem_error],
    [:CA1,:PFC,:CA1PFC])

    @info "loop" datacut region

    datacutStr = string(datacut)

    spikes, beh = Load.load("RY16",36; data_source=["spikes","behavior"])
    @time i = I[bestpartialmatch(keys(F), (;datacut, widths=5))];
    @time f = F[bestpartialmatch(keys(F), (;datacut, widths=5))];
    f = ShiftedFields(deepcopy(f))
    unitshift = Timeshift.types.matrixform(f)
    getshift(arrayOfFields::DimArray, s) = arrayOfFields[:, arrayOfFields.dims[2].==s];
    shifts = collect(unitshift.dims[2])
    push_celltable!( unitshift, cells, :unit, :area)
    push_dims!(unitshift)
    push_shiftmetric!(unitshift, best_tau!; metric=:bitsperspike)
    #push_dims!(f, vec(getindex.(m, :coherence)); dim=:unit, metric=:coh_at_zero)

    # Which cells pass our criteria?
    metricfilter = metricfilters[region]

    # Get filtered shift=0 fields
    shift0 = filter(metricfilter, getshift(unitshift,0))
    shift0 = shift0[sortperm(shift0[:unit])]

    spikes = subset(spikes, :unit=>x->Utils.squeeze(any(x .∈ shift0[:unit]',dims=2)))
    spikes = Utils.filtreg.filterAndRegister(beh,spikes; filters=filt[datacut], on="time",  transfer=["velVec"], filter_skipmissingcols=true)[2]
    sort(unique(spikes.unit))

    # ------------------------
    # Firing rate correlations
    # ------------------------

    # Snage a rate matrix (time × neuron)
    R = torate(spikes , beh; gaussian=0)
    # Zscore
    R = (R .- Array(mean(R, dims=:time)))./Array(std(R,dims=:time))

    heatmap(R', clim=(-3, 3))
    shift0_filtbyrate = shift0[vec(any(shift0[:unit] .∈ collect(R.dims[2])', dims=2))]

    # Compute correlation
    Rcor = cor(R, dims=:time)
    #heatmap(Array(Rcor), c=:vik, clim=(-1,1))

    # ------------------------
    # Unit τ differences Δ
    # ------------------------

    τ = NaN .* ones(size(Rcor))
    (u1,u2,r) = first(zip(Rcor.dims[1], Rcor.dims[2], Rcor))
    for (i1,i2) in Iterators.product(1:size(Rcor,1), 1:size(Rcor,2))
        u1, u2 = Rcor.dims[1].val[i1],
                 Rcor.dims[2].val[i2]
        if u1 == u2
            continue
        end
        τ1, τ2 = shift0_filtbyrate[findfirst(shift0_filtbyrate[:unit] .== u1)][:bestshift_bitsperspike],
                 shift0_filtbyrate[findfirst(shift0_filtbyrate[:unit] .== u2)][:bestshift_bitsperspike]
        τ[i1,i2] = τ1 - τ2
    end


    # ------------------------
    # Rate corr ∝ Δτ
    # ------------------------
    vecTau = vec(abs.(τ))
    vecCor = vec(abs.(Rcor))
    good = (!).(isnan.(vecTau))
    vecTau, vecCor = vecTau[good], vecCor[good]
    df = DataFrame(Dict(:vecTau=>vecTau, :vecCor=>vecCor))
    df.rank = sortperm(df.vecCor)
    df.bintau = Int.(floor.(Utils.norm_extrema(df.vecTau, [0,4])))
    combine(groupby(df, :bintau),nrow,:vecTau=>mean)
    #transform!(df, :vecCor => (x->log10.(x)) => :vecCor)

    Plot.setfolder("timeshift", "functional_connectivity")

    q = 95
    linmod = fit(LinearModel, @formula(vecCor ~ vecTau), df)
    fval=round(ftest(linmod.model).pval; sigdigits=2)
    r2stat = round(r2(linmod); sigdigits=2)
    @df df scatter(:bintau, :vecCor, label="samples", xlabel="δτ=τ₁-τ₂", ylabel="pearson", title="Datacut=$datacut\nr2=$r2stat, pval=$fval")
    dfs = combine(groupby(df, :bintau), :vecCor=>(x->quantile(x,q/100)))
    @df dfs plot!(:bintau, :vecCor_function, label="Q_$q")
    dfs = combine(groupby(df, :bintau), :vecCor=>mean)
    @df dfs plot!(:bintau, :vecCor_mean, label="mean")
    ŷ = predict(linmod, df)
    plot!(df.vecTau, ŷ, c=:black)
    Plot.save((;desc="linear model", datacut=datacutStr, region))
    linmod

    q = 95
    linmod = fit(LinearModel, @formula(rank ~ vecTau), df)
    fval=round(ftest(linmod.model).pval; sigdigits=2)
    r2stat = round(r2(linmod); sigdigits=2)
    @df df scatter(:bintau, :rank, label="samples", xlabel="δτ=τ₁-τ₂", ylabel="pearson", title="Datacut=$datacut\nr2=$r2stat, pval=$fval")
    dfs = combine(groupby(df, :bintau), :rank=>(x->quantile(x,q/100)))
    @df dfs plot!(:bintau, :rank_function, label="Q_$q")
    dfs = combine(groupby(df, :bintau), :rank=>mean)
    @df dfs plot!(:bintau, :rank_mean, label="mean")
    ŷ = predict(linmod, df)
    plot!(df.vecTau, ŷ, c=:black)
    Plot.save((;desc="rank linear model", datacut=datacutStr, region))
    linmod

    #alpha=0.4
    #ecdfplot( @subset(df, :bintau .== 0).vecCor; fillrange=0, alpha, label="1st bin")
    #ecdfplot!(@subset(df, :bintau .== 1).vecCor; fillrange=0, alpha, label="2nd bin")
    #ecdfplot!(@subset(df, :bintau .== 2).vecCor; fillrange=0, alpha, label="3rd bin", legend=:outerbottomright)
    #ecdfplot!(@subset(df, :bintau .== 3).vecCor; alpha, fillrange=0, label="4th bin", legend=:outerbottomright)
    #ecdfplot!(@subset(df, :bintau .== 4).vecCor; alpha, fillrange=0, label="5th bin", legend=:outerbottomright)

    #@df df violin(:vecCor, group=:bintau)
    #
    #import Gadfly
    #p = Gadfly.plot(x=vecTau, y=vecCor,
    #            Gadfly.Stat.binmean, Gadfly.Geom.point, Gadfly.Geom.line)


    # ---------------
    # XY CORRELATIONS
    # ---------------
    Plot.setfolder("timeshift", "xy", "centroids")

    c = get.([ColorSchemes.vik],
            Utils.norm_extrema(shift0[:bestshift_bitsperspike], [0,1]))
    task.plotboundary(taskdf, seriestype=:shape, alpha=0.2)
    scatter!(eachcol(hcat(shift0[:centroid]...)')...; c, label="", framestyle=:none, title="Field centroids annotated betst τ")
    Plot.save((;desc="Best_tau_of_field_centroids", datacut=datacutStr, region))

       

    # -------------------------
    # FUTURE/PRESENT/PAST INFO
    # shift = 0
    # -------------------------
    Plot.setfolder("timeshift","xy","multiunit_quintiles")

    # Create quintiles of bestshift and plot field averages
    shift0[:bestshift_bins] =  Utils.binning.digitize(shift0[:bestshift_bitsperspike], 5)

    averages_t0 = Dict(
        :standard=>[],
        :norm01 =>[],
        :normBSP=>[],
        :normB=>[],
        :N=>[],
        )
    for bin = 1:5
        fields_in_bin = shift0[shift0[:bestshift_bins] .== bin]
        to_average = [f.rate for f in fields_in_bin]
        bitsperspike  = fields_in_bin[:bitsperspike]
        bits = fields_in_bin[:bitsperspike] .* fields_in_bin[:totalcount]
        if isempty(to_average)
            @warn "empty to_average" datacut region bin
        end
        try
            push!(averages_t0[:standard], mean(to_average))
        catch
            @infiltrate
        end
        push!(averages_t0[:norm01],   mean(map(x->Utils.nannorm_extrema(x, [0, 1]),
                                            to_average)))
        bitsperspike_norm(x::AbstractArray,bsp::Real) = Utils.nannorm_extrema(x, [0, bsp])
        push!(averages_t0[:normBSP],  mean(map(bitsperspike_norm,
                                            to_average, bitsperspike)))
        push!(averages_t0[:normB],  mean(map(bitsperspike_norm,
                                            to_average, bits)))
        push!(averages_t0[:N], length(to_average))
    end

    # Field averages, without scaling each neuron
    plot([heatmap(aver, aspect_ratio=1, framestyle=:none,
                  title="Best, standard\n$datacut, N=$(averages_t0[:N])") 
          for aver in averages_t0[:standard]]...)
    Plot.save((;datacut=datacutStr, shift=0, average="standard", region))

    # Field averages, scaling each neuron 0-1
    plot([heatmap(aver, aspect_ratio=1, framestyle=:none,
                  title="Best, standard\n$datacut, N=$(averages_t0[:N])") 
          for aver in averages_t0[:norm01]]...)
    Plot.save((;datacut=datacutStr, shift=0, average="norm01", region))

    # Field averages, scaling each neuron by bits per spike
    plot([heatmap(aver, aspect_ratio=1, framestyle=:none,
                  title="Best, standard\n$datacut, N=$(averages_t0[:N])") 
          for aver in averages_t0[:normBSP]]...)
    Plot.save((;datacut=datacutStr, shift=0, average="normBSP", region))

    plot([heatmap(aver, aspect_ratio=1, framestyle=:none,
                  title="Best, standard\n$datacut, N=$(averages_t0[:N])") 
          for aver in averages_t0[:normB]]...)
    Plot.save((;datacut=datacutStr, shift=0, average="normB", region))

    # -------------------------
    # FUTURE/PRESENT/PAST INFO
    # shift = best
    # -------------------------
    Plot.setfolder("timeshift","xy","multiunit_quintiles")

    # Create quintiles of bestshift and plot field averages
    B = [findfirst(row) for row in eachrow(unitshift[:,1][:bestshift_bitsperspike].data .== shifts')]
    shiftB = [unit[b] for (unit,b) in zip(eachrow(unitshift),B)]
    shiftB = filter(metricfilter, shiftB)

    shiftB[:bestshift_bins] =  Utils.binning.digitize(
                                shiftB[:bestshift_bitsperspike], 5)

    averages_tB = Dict(
        :standard=>[],
        :norm01 =>[],
        :normBSP=>[],
        :normB=>[],
        :N=>[],
        )
    for bin = 1:5
        fields_in_bin = shiftB[shiftB[:bestshift_bins] .== bin]
        to_average = [f.rate for f in fields_in_bin]
        if isempty(to_average)
            @warn "empty to_average" datacut region bin
        end
        bitsperspike  = fields_in_bin[:bitsperspike]
        bits = fields_in_bin[:bitsperspike] .* fields_in_bin[:totalcount]
        push!(averages_tB[:standard], mean(to_average))
        push!(averages_tB[:norm01],   mean(map(x->Utils.nannorm_extrema(x, [0, 1]),
                                            to_average)))
        bitsperspike_norm(x::AbstractArray,bsp::Real) = Utils.nannorm_extrema(x, [0, bsp])
        push!(averages_tB[:normBSP],  mean(map(bitsperspike_norm,
                                            to_average, bitsperspike)))
        push!(averages_tB[:normB],  mean(map(bitsperspike_norm,
                                            to_average, bits)))
        push!(averages_tB[:N], length(to_average))
    end

    # Field averages, without scaling each neuron
    plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
          for aver in averages_tB[:standard]]..., title="Best, standard\n$datacut, N=$(averages_tB[:N])")
    Plot.save((;datacut=datacutStr, shift="best", average="standard", region))

    # Field averages, scaling each neuron 0-1
    plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
          for aver in averages_tB[:norm01]]..., title="Best, norm01\n$datacut, N=$(averages_tB[:N])")
    Plot.save((;datacut=datacutStr, shift="best", average="norm01", region))

    # Field averages, scaling each neuron by bits per spike
    plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
          for aver in averages_tB[:normBSP]]..., title="Best, normBSP\n$datacut, N=$(averages_tB[:N])")
    Plot.save((;datacut=datacutStr, shift="best", average="normBSP", region))

    # Field averages, scaling each neuron by bits per spike
    plot([heatmap(aver, aspect_ratio=1, framestyle=:none) 
          for aver in averages_tB[:normB]]..., title="Best, normB\n$datacut, N=$(averages_tB[:N])")
    Plot.save((;datacut=datacutStr, shift="best", average="normB", region))

    # -------------------------
    # FUTURE/PRESENT/PAST INFO
    #
    # Rather than *just* best shifts, 
    # we can also try all shifts
    # -------------------------
    Plot.setfolder("timeshift","xy","multiunit_allCellsVote")

    multiunit = Vector{AbstractArray}(undef, size(unitshift,2))
    for (s,shift) in enumerate(eachcol(unitshift))
        # Get multiunit response weighted by bits
        shift = filter(metricfilter, shift)
        bitsperspike, totalcount = shift[:bitsperspike],
                                   shift[:totalcount]
        bits = bitsperspike .* totalcount
        units = [x.rate for x in shift]
        multiunit[s] = mean(units .* bits)
    end
    c = get(ColorSchemes.vik, Utils.norm_extrema(shifts,[0,1]))
    mx = maximum(nanmaximum.(multiunit))
    mi = minimum(nanminimum.(multiunit))
    anim = @animate for s in 1:length(multiunit)
        heatmap(multiunit[s], aspect_ratio=1, title="population, bits weighted, $datacut\n$(unitshift[1,s][:shift])", titlefonthalign=:center, titlefontcolor=c[s], clim=(mi,mx))
    end
    desc = replace(string((;datacut)),"("=>"", ")" => "", "\""=>"", ","=>"", " " => "", ":" => "")
    filename = plotsdir("timeshift","xy","multiunit_allCellsVote", "$(desc)_$region.gif")
    gif(anim, filename; fps=10)

end
