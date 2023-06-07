# -------
# INFO:
# -------
# - iso_* : a dataframe with summarized information
# - per_* : a dataframe with information about each period
# --------
# BUG:
# --------
# DO NOT trust any non-`period` field containing call to get_isolation_summary
# just yet!
# ---------
# TODO:
# Shantanu
# ---------
# - [ ] Split by cuemem
# - [ ] Split by interneuron
# - [ ] Split by ha

opt = 
    Dict("animal"       => "RY16",
        "day"           => 36,
        "checkpoint"    => true,
        "tetrode_set"   => "synced",
        "filt"          => :all
    )
tetrode_sets = ["synced", "best"]

for ((animal,day), tetrode_set) in Iterators.product(DI.animaldays(),
                                                    tetrode_sets)

    @eval Main include(scriptsdir("isolated", "load_isolated.jl"))
    # @eval Main include("./load_isolated.jl")
    @assert :isolated ∈ propertynames(spikes)
    @assert beh !== nothing
    Munge.nonlocal.setunfilteredbeh(beh)

    using DI.Labels, DataFrames, DataFramesMeta, StatsPlots
    using DIutils.plotutils
    using GoalFetchAnalysis.Munge.isolated
    using GoalFetchAnalysis.Munge.nonlocal
    using ElectronDisplay

    Plot.setparentfolder("isolated")
    Plot.setappend("$animal-$day"-"$tetrode_set")

    # Add a column to our spikes dataframe about its cell's meanrate
    DIutils.filtreg.register(cells, spikes, on="unit", transfer=["meanrate", "area"])
    # Add a behavioral info to spikes
    DIutils.filtreg.register(beh,   spikes, on="time", transfer=["velVec", 
        "period","correct", "cuemem",
        "hatraj", "ha"])
    @assert :period ∈ propertynames(spikes)

    # Acquire the isolation table
    iso_sum  = get_isolation_summary(spikes)

    # Obtain the isolated spikes dataframe
    isol = last(groupby(subset(spikes, :isolated=>x->(!).(isnan.(x)), skipmissing=true),
                                :isolated))
    @assert all(isol.isolated .== true)

    # Select interneuron, roughly
    spikes.interneuron = spikes.meanrate .> 5
    histogram(cells.meanrate)
    #   BUG: MISSING NOT FOUND 👈︎
    iso_ct = get_isolation_summary(spikes,[:ha, :cuemem, :interneuron])
    sort!(iso_ct, [:area, :interneuron, :cuemem])

    per_ct = get_isolation_summary(spikes, [:ha, :cuemem, :interneuron, :period, :correct])
    @subset!(per_ct, :events_per_time .!= Inf)
    sort!(per_ct, [:area, :interneuron, :cuemem, :period])
    @subset!(per_ct, (:cuemem .== -1 .&& :correct .== -1) .||
                                   (:cuemem .== 0 .&& :correct .!= -1)  .||
                                   (:cuemem .== 1 .&& :correct .!= -1))
    subset!(per_ct, :events_per_time => x-> (!).(isinf.(x)))
    per_ct


    # ==============================
    # SECITON : BASIC ISO SPIKING
    # ==============================
    # Just the concentration of spikes at phase with and without iso
    DIutils.filtreg.register(lfp, spikes, on="time", transfer=["phase"])

    # ALL SPIKES
    sp = dropmissing(spikes, :isolated)
    p1 = @df @subset(sp, :isolated .== false) histogram(:phase, bins=range(0,stop=2pi,length=20), normed=true, 
        label=["non-isolated"], legend=:topleft)
    p2 = @df @subset(sp, :isolated .== true) histogram(:phase, bins=range(0,stop=2pi,length=20), normed=true, 
        label=["isolated"], legend=:topleft, xlabel="Phase")
    plot(p1,p2,layout=(2,1), size=(800,600), legend=:topleft, legendfontsize=10, 
        legendtitle="Status", title="$animal-$day-$tetrode_set")
    # INDIVIDUAL CELLS
    using DirectionalStatistics
    using IntervalSets
    spc = combine(groupby(dropmissing(spikes,:isolated), [:unit,:isolated]),
        :phase => (x->Circular.mean(x,0..2pi)) => :mean_phase,
    )
    p1 = @df @subset(spc, :isolated .== false) histogram(:mean_phase, bins=30, normed=true, 
        label=["non-isolated"], legend=:topleft)
    p2 = @df @subset(spc, :isolated .== true) histogram(:mean_phase, bins=30, normed=true, 
        label=["isolated"], legend=:topleft, xlabel="Phase")
    plot(p1,p2,layout=(2,1), size=(800,600), legend=:topleft, legendfontsize=10, legendtitle="Status", title="$animal-$day-$tetrode_set")

end


#     _  _     _____ ___  ____   ___  ____  
#  _| || |_  |_   _/ _ \|  _ \ / _ \/ ___| 
# |_  ..  _|   | || | | | | | | | | \___ \ 
# |_      _|   | || |_| | |_| | |_| |___) |
#   |_||_|     |_| \___/|____/ \___/|____/
#                                          
# - labels
#     - titles
#     - legend titles
#     - any missing x/ylabels
# - more accurate peak to peak
# - split these by moving and still
# - split by 1st and 2nd visit
#

# ==============================
# SECTION : ALL SPIKES
# NOTE: Not period wise, so less valuable, (•___•)
#    _  _        _    _     _       ____  ____ ___ _  _______ ____  
#  _| || |_     / \  | |   | |     / ___||  _ \_ _| |/ / ____/ ___| 
# |_  ..  _|   / _ \ | |   | |     \___ \| |_) | || ' /|  _| \___ \ 
# |_      _|  / ___ \| |___| |___   ___) |  __/| || . \| |___ ___) |
#   |_||_|   /_/   \_\_____|_____| |____/|_|  |___|_|\_\_____|____/ 
# SHEER DIFFFERENT RATE 📈 OF EVENTS
# (adjacent/isolated) over cue, mem, nontask
# ==============================
# =========PLOTS === NOT SPLIT BY CELL TYPE ============
#begin
Plot.setfolder("MUA and isolated MUA", "allspikes")

    Plot.setfolder("MUA and isolated MUA")
    kws=(;legend_position=Symbol("outerbottomright"))
    @df @subset(iso_sum,:area.=="CA1") bar(:cmlab, :timespent, ylabel="Time spent", group=:cuemem; kws...,
    legendtitle="Task", fillalpha=0.5)
    Plot.save((;desc="time spent"))
    @df iso_sum bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:cuemem; kws..., 
                   legendtitle="Task", fillalpha=0.5)
    Plot.save((;desc="MUA per second"))
    @df iso_sum bar(:cuearea, :isolated_mean, ylabel="Isolated Spikes (sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:cmlab; kws..., 
                   legendtitle="Task")
    Plot.save((;desc="fraction of isolated spikes"))
    @df iso_sum bar(:cuearea, :isolated_events_per_time, ylabel="Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:cuemem; kws..., 
                   legendtitle="Task")
    Plot.save((;desc="isolated spikes per second"))

#end
# ======================================================

# =========PLOTS ===SPLIT BY CELL STATE 🔺🔺/ 🔵 =======

# ISSUE: HAVE TO MAKE SURE THAT NONTASK CUE AND MEM TIMESPENT ARE COUNTED CORRECTLY

unique_interneurons = unique(iso_ct.interneuron)
yvals = [(:events_per_time, "MUA events per second\n$(filt_desc[:all])"),
         (:isolated_mean, "Isolated Spikes / All spikes\nsign of CA1-PFC interaction)\n$(filt_desc[:all])"),
         (:isolated_events_per_time, "Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])")]
P = OrderedDict()
interneuron = unique_interneurons[1]
for interneuron in unique_interneurons
    df_subset = filter(row -> row.interneuron == interneuron, iso_ct)
    (yval, ylabel) = yvals[1]
    for (yval, ylabel) in yvals
        p = bar(df_subset.areacuememha, df_subset[!,yval]; ylabel=ylabel, alpha=0.5, size=(1000,500), 
            leftmargin=5*Plots.mm, group=df_subset.ha, bottommargin=5*Plots.mm, title=neuron_names[interneuron])
        # Replace `savefig` with the correct function to save your plot
        Plot.save("celltype_$(interneuron)_$(yval).png")
        P[interneuron, yval] = p
    end
end
plot(values(P)..., size=(1000,500).*1.2)
Plot.save("celltype_all.png")



# ======================================================
# NOTE:  🤑🤑🤑
#    _  _     ____           _           _               _          
#  _| || |_  |  _ \ ___ _ __(_) ___   __| |    __      _(_)___  ___ 
# |_  ..  _| | |_) / _ \ '__| |/ _ \ / _` |____\ \ /\ / / / __|/ _ \
# |_      _| |  __/  __/ |  | | (_) | (_| |_____\ V  V /| \__ \  __/
#   |_||_|   |_|   \___|_|  |_|\___/ \__,_|      \_/\_/ |_|___/\___|
#                                                                   

# =========PLOTS =======================================
# -------- mua per second ------------------------------

# CELL AGNOSTIC ⬜

    # -------- percent isolated spikes ---------------------
    # color by cue/mem
    Plot.setfolder("MUA and isolated MUA", "period-wise_cuemem_pyrint")
    unique_interneurons = unique(per_ct.interneuron)
    unique_areas = unique(iso_ct.area)
    yvals = [(:events_per_time, "mua events per second\n$(filt_desc[:all])"),
             (:isolated_mean, "isolated spikes / all spikes\nsign of ca1-pfc interaction)\n$(filt_desc[:all])"),
             (:isolated_events_per_time, "isolated mua × sec⁻1\n(sign of ca1-pfc interaction)\n$(filt_desc[:all])")]
    P = OrderedDict()
    for (interneuron, area) in Iterators.product(unique_interneurons, unique_areas)
        for yval in yvals
            df_subset = filter(row -> row.interneuron == interneuron && row.area == area, per_ct)
            # quantile filter yval
            if isempty(df_subset)
                continue
            end
            df_subset = subset(df_subset, yval[1]=>row->row .< quantile(row|>DIutils.skipnan, 0.98))
            if isempty(df_subset)
                continue
            end
            p = scatter(df_subset.cuemem .+ 0.25.*(df_subset.correct .- 0.5) .+ (randn(size(df_subset.cuemem)) .* 0.05),
            df_subset[!,yval[1]], group=df_subset.correct, alpha=0.4,
                        xtick=([-1,0,1],["nontask","cue","mem"]),
                        title="$(neuron_names[interneuron]) $(area)",
                        xlabel="task", ylabel=yval[2], 
                        legend_title="correct/error", legend_position=:outerbottomright)
            boxplot!(df_subset.cuemem .+ 0.25.*(df_subset.correct .- 0.5),
            df_subset[!,yval[1]], alpha=0.5, linewidth=2, linecolor=:black,
            group=df_subset.correct, legend=false, outlier=false)
            # replace `savefig` with the correct function to save your plot
            Plot.save("isolation_$(interneuron)_$(yval[1])_$(area).png")
            P[(interneuron, area, yval[1])] = p
        end
    end
    plot(values(P)..., size=(1000,500).*1.4, leftmargin=5*Plots.mm, bottommargin=5*Plots.mm)
    Plot.save("metrics_cuemem_pyrint_all.png")
    current()

    Plot.setfolder("MUA and isolated MUA", "period-wise_cuememha_pyrint")
    unique_interneurons = unique(per_ct.interneuron)
    unique_areas = unique(iso_ct.area)
    yvals = [(:events_per_time, "mua events per second\n$(filt_desc[:all])"),
     (:isolated_mean, "isolated spikes / all spikes\nsign of ca1-pfc interaction)\n$(filt_desc[:all])"),
     (:isolated_events_per_time, "isolated mua × sec⁻1\n(sign of ca1-pfc interaction)\n$(filt_desc[:all])")]
    P = OrderedDict()
    (interneuron, area) = (unique_interneurons[1], unique_areas[1])
    for (interneuron, area) in Iterators.product(unique_interneurons, unique_areas)
        yval = yvals[1]
        for yval in yvals
            df_subset = filter(row -> row.interneuron == interneuron && row.area == area, per_ct)
            df_subset = subset(df_subset, yval[1] => row-> row .< quantile(row|>DIutils.skipnan, 0.99))
            if isempty(df_subset)
                continue
            end
            p = scatter(df_subset.cuememhaind .+ 0.25.*(df_subset.correct .- 0.5) .+ (randn(size(df_subset.cuememhaind)) .* 0.05),
                        df_subset[!,yval[1]], group=df_subset.correct, alpha=0.4,
                        xtick=([1,2,3,4,5], ["AC", "HC", "HM", "AM", "A-"]),
                        title="$(neuron_names[interneuron]) $(area)",
                        xlabel="task", ylabel=yval[2],
                        label="")
            # boxplot!(df_subset.cuememhaind .+ 0.25.*(df_subset.correct .- 0.5),
            #          df_subset[!,yval[1]], alpha=0.5, linewidth=2, linecolor=:black,
            #          group=df_subset.correct, legend=false, outlier=false)
            # Replace `savefig` with the correct function to save your plot
            Plot.save("metrics_$(interneuron)_$(yval[1])_$(area)")
            P[(interneuron, area, yval[1])] = p
        end
    end
    plot(values(P)..., size=(1000,500).*1.4, leftmargin=5*Plots.mm, bottommargin=5*Plots.mm)
    Plot.save("metrics_cuememha_pyrint_all")


    for (interneuron, area) in Iterators.product(unique_interneurons, unique_areas)
        for yval in yvals
            df_subset = filter(row -> row.interneuron == interneuron && row.area == area, per_ct)
            df_subset = subset(df_subset, yval[1]=>row->row .< quantile(row|>DIutils.skipnan, 0.99))
            if isempty(df_subset)
                continue
            end
            p = P[(interneuron, area, yval[1])]
            plot(p)
            boxplot!(df_subset.cuememhaind .+ 0.25.*(df_subset.correct .- 0.5),
                     df_subset[!,yval[1]], alpha=0.5, linewidth=2, linecolor=:black,
                     group=df_subset.correct, legend=false, outlier=false)
            # Replace `savefig` with the correct function to save your plot
            Plot.save("boxplot_metrics_$(interneuron)_$(yval[1])_$(area)")
            P[(interneuron, area, yval[1])] = p
        end
    end
    plot(values(P)..., size=(1000,500).*1.4, leftmargin=5*Plots.mm, bottommargin=5*Plots.mm)
    Plot.save("boxplot_metrics_cuememha_pyrint_all")


# Animations 🎥
    # -------- multi - metric : one predict other ? --------
    Plot.setfolder("MUA and isolated MUA", "period-wise_cuemem_pyrint")
    for link_axes in [false, true]
        anim = @animate for i in 1:360
            p1=@df @subset(per_ct, :interneuron .== false, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                        group=:correct, alpha=0.5, camera=(i,30), 
                        title="pyr cue",
                        xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s")
            p2=@df @subset(per_ct, :interneuron .== false, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                        group=:correct, alpha=0.5, camera=(i,30),
                        title="pyr mem",
                        xlabel="fraction(isolated)", ylabel="isolated multiunit/s", zlabel="multiunit/s")
            plot(p1,p2; link = link_axes ? :all : :none)
        end
        gif(anim, plotsdir(Plot.complete_folder_args..., "separatrix-three-vars-gif_pyr_facet=cuemem_group=correct_link=$link_axes.gif"); loop=1)
        anim = @animate for i in 1:360
            p1=@df @subset(per_ct, :interneuron .== false, :cuemem .== 1) begin
                scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr mem")
            end
            p2=@df @subset(per_ct, :interneuron .== true, :cuemem .== 1) begin
                scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30),ylabel="isolated multiunit/s", zlabel="multiunit/s", title="int mem")
            end
            p3=@df @subset(per_ct, :interneuron .== false, :cuemem .== 0) begin
                scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr cue")
            end
            p4=@df @subset(per_ct, :interneuron .== true, :cuemem .== 0) begin
                scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30),ylabel="isolated multiunit/s", zlabel="multiunit/s", title="int cue")
            end
            plot(p3,p4,p1,p2; layout=grid(2,2), link = link_axes ? :all : :none)
        end
        gif(anim, plotsdir(Plot.complete_folder_args..., "separatrix-three-vars-gif_facet=pyrint,cuemem_group=correct_link=$link_axes.gif"); loop=1)
    end

# ======================================================

Plot.setfolder("MUA and isolated MUA", "period-wise-summary")
#=
===========================
summary of period iso_mean
===========================
=#


P=[]
# ISSUE: FIX 
plot_settings = [(false, "CA1", "CA1 pyr"), (true, "CA1", "CA1 int"), (false, "PFC", "PFC pyr")]
for (interneuron_value, area_value, title_value) in plot_settings
    XX = combine(groupby(@subset(per_ct, :interneuron .== interneuron_value, :area .== area_value),
                  :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x)))
    kws=(;xlabel="cuemem", ylabel="isolated frac", title=title_value)
    @df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function; kws...)
    Plot.save(kws)
    push!(P, p)
end
plot(P..., size=(1000,500).*1.4, leftmargin=5*Plots.mm, bottommargin=5*Plots.mm)
Plot.save("isomean_summary_cuememha_pyrint_all")


#================================
summary of period mua
=================================#
P = []
plot_settings = [(nothing, "CA1", "CA1 MUA per second"), (false, "CA1", "CA1 pyr MUA"), 
         (true, "CA1", "CA1 int MUA"), (false, "PFC", "PFC pyr")]
for (interneuron_value, area_value, title_value) in plot_settings
    if interneuron_value === nothing
        XXg = groupby(@subset(per_ct, :area .== area_value), :cuemem)
    else
        XXg = groupby(@subset(per_ct, :interneuron .== interneuron_value, :area .== area_value), :cuemem)
    end
    XX = combine(XXg, :events_per_time => median, 
                 :events_per_time => x->nanstd(x)./sqrt(length(x)))
    kws=(;xlabel="cuemem", ylabel="MUA per sec", title=title_value)
    @df XX bar(:cuemem, :events_per_time_median, yerr=:events_per_time_function; linewidth=2, kws...)
    @df XX plot!(:cuemem, :events_per_time_median, yerr=:events_per_time_function, linestyle=:none)
    Plot.save(kws)
end
plot(P..., size=(1000,500).*1.4, leftmargin=5*Plots.mm, bottommargin=5*Plots.mm)
Plot.save("mua_summary_cuememha_pyrint_all")


#============================
summary of iso events per time
============================#

P = []
parameters = [(false, "CA1", "CA1 pyr"), (true, "CA1", "CA1 int"), (false, "PFC", "PFC pyr")]
for (interneuron_value, area_value, title_value) in parameters
    XX = combine(
        groupby(@subset(per_ct, :interneuron .== interneuron_value, :area .== area_value), :cuemem), 
        :isolated_events_per_time => median, 
        :isolated_events_per_time => x->std(x)/sqrt(length(x))
    )
    kws=(;xlabel="cuemem", ylabel="iso events per sec", title=title_value)
    @df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function; kws...)
    Plot.save(kws)
end
plot(P..., size=(1000,500).*1.4, leftmargin=5*Plots.mm, bottommargin=5*Plots.mm)
Plot.save("isoevents_summary_cuememha_pyrint_all")

"""
Mean cycle dist
"""

Plot.setfolder("isolation-basic_column_plots")
histogram(filter(x->!ismissing(x) && x<100, @subset(spikes,:area.=="CA1", :interneuron .!= true).meancyc), bins=40, yscale=:log10)
Plot.save((;desc="Mean cycle, all cells"))
histogram(filter(x->!ismissing(x) && x<100, @subset(spikes,:area.=="CA1", :interneuron .!= true).nearestcyc), bins=40, yscale=:log10)
Plot.save((;desc="Nearest cycle, all cells"))

"""
===========================
SVM: All three variables
===========================
"""
Plot.setfolder( "svm")

    using MLJ
    MLJ.@load SVC pkg=LIBSVM
    cv = StratifiedCV(nfolds=2, shuffle=true)
    #cv = Holdout(fraction_train=0.80, shuffle=true)

    function eva(cv::ResamplingStrategy, svc::Machine; n=100)
        measure=[accuracy, measure_0, measure_1]
        #if typeof(cv) <: StratifiedCV
        #    res = MLJ.evaluate!(svc; resampling=cv,
        #                       verbosity=0,measure)
        #    df = DataFrame(res)
        #elseif typeof(cv) <: Holdout
            res = [MLJ.evaluate!(svc; resampling=cv,
                               verbosity=0, measure)
                  for _ in 1:n]
            df = [transform!(DataFrame(r), :measure=>(x->i*ones(size(DataFrame(r),1)))=>:replicate) 
                  for (i,r) in enumerate(res)]
            df = vcat(df...)
        #end
        res, df
    end
    function measure_0(y,ŷ)
        inds = y .== 0
        accuracy(y[inds], ŷ[inds])
    end
    function measure_1(y,ŷ)
        inds = y .== 1
        accuracy(y[inds], ŷ[inds])
    end
    function get_confusion_matrix()
    end
    cmt=Dict("accuracy"=>"accuracy", "measure_0 "=>"pred. cue\npred. error", "measure_1 "=>"pred. memory\npred. correct")

    svm_data = @subset(per_ct, :correct .== 0 .|| :correct .== 1);
    svm_data = dropmissing(svm_data[:,[:isolated_mean,:events_per_time,:isolated_events_per_time,:correct, :cuemem]]);
    XX = Matrix(svm_data[!, [:isolated_mean,:events_per_time,:isolated_events_per_time]]);

    # ------------------------------------------------------------
    # KEY QUESTION: Can we predict CORRECT/ERROR using ALL 3 
    # [:isolated_mean,:events_per_time,:isolated_events_per_time]
    # ------------------------------------------------------------
        y = categorical(svm_data[!,:correct]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        ld = DIutils.mlj.measureidxDict(df)
        xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
        @df df  begin
            scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)),
                    :measurement; ylim=(0,1), label="Predict correct", xtick)
        end

    # ------------------------------------------------------------
    # KEY QUESTION: Can we predict CUE-MEM using ALL 3 
    # [:isolated_mean,:events_per_time,:isolated_events_per_time]
    # ------------------------------------------------------------
        y = categorical(svm_data[!,:cuemem]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
        @df df begin
            scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)),
                     :measurement; ylim=(0,1), label="Predict cuemem", legend=:none, xtick)
        end
    Plot.save("all three vars")

    """
    ===========================
    SVM: events_per_time
    ===========================
    """

        svm_data = @subset(per_ct, :correct .== 0 .|| :correct .== 1);
        svm_data = dropmissing(svm_data[:,[:events_per_time,:correct, :cuemem]]);
        XX = MLJ.table(Matrix(svm_data[!, [:events_per_time,]]))

        y = categorical(svm_data[!,:correct]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        xtick = (collect(values(ld)), getindex.([cmt], collect(keys(ld))))
        @df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
                       xtick)

        y = categorical(svm_data[!,:cuemem]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        @df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
                       xtick)

    Plot.save("events per time")

    """
    ===========================
    SVM: isolated_mean
    ===========================
    """

        svm_data = @subset(per_ct, :correct .== 0 .|| :correct .== 1);
        svm_data = dropmissing(svm_data[:,[:isolated_mean,:correct, :cuemem]]);
        XX = MLJ.table(Matrix(svm_data[!, [:isolated_mean,]]))

        y = categorical(svm_data[!,:correct]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
        @df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict correct",
                       xtick)

        y = categorical(svm_data[!,:cuemem]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
        @df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)), :measurement; ylim=(0,1), label="Predict cuemem",
                       xtick)

        Plot.save("isolation")

    """
    ===========================
    SVM: isolated_events_per_time
    ===========================

    """

        svm_data = @subset(per_ct, :correct .== 0 .|| :correct .== 1);
        svm_data = dropmissing(svm_data[:,[:isolated_events_per_time,:correct, :cuemem]]);
        XX = MLJ.table(Matrix(svm_data[!, [:isolated_events_per_time,]]))

        y = categorical(svm_data[!,:correct]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
        @df df scatter(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)),
                       :measurement; ylim=(0,1), label="Predict correct", xtick)

        y = categorical(svm_data[!,:cuemem]);
        svc_mdl = MLJLIBSVMInterface.SVC(degree=Int32(4));
        svc = machine(svc_mdl, XX, y);
        res, df = eva(cv, svc)
        xtick = (collect(values(ld)), getindex.([cmt],collect(keys(ld))))
        @df df scatter!(:idxmeasure .+ 0.1 .* randn(size(:idxmeasure)),
                        :measurement; ylim=(0,1), label="Predict cuemem", xtick)
        Plot.save("isolated_events_per_time")


"""
===========================
Relation of isolated to out of field?
===========================

Answer: Confusing or not much! At least explored with hulls. Maybe still fits Jai's
"""

isoin_sum = combine(groupby(spikes, [:isolated,:interneuron]),:infield => x->mean(skipmissing(x)), renamecols=false)
dropmissing!(isoin_sum)

@df isoin_sum bar(:isolated, :infield, group=:interneuron; xticks=([0,1],["adj","iso"]), legend_title="IsInterneuron")

isoin_sum_celltype = combine(groupby(spikes, [:unit, :isolated, :interneuron]),:infield=>x->nanmean(skipmissing(x));renamecols=false)
isoin_sum_celltype[!,:infield] = replace(isoin_sum_celltype.infield, NaN=>missing)
dropmissing!(isoin_sum_celltype)

Plot.setfolder("isolation-infield")
@df @subset(isoin_sum_celltype,:interneuron .!= true) boxplot(:isolated, :infield, xticks=([0,1],["adj","iso"]))
@df @subset(isoin_sum_celltype,:interneuron .!= true) scatter!(:isolated + randn(size(:isolated)).*0.1, :infield, xticks=([0,1],["adj","iso"]))
Plot.save((;desc="Not much difference adjacent-iso in percent infield spikes",hullmax=1,thresh=0.8))


"""
===========================
Isolated spikes ISI
===========================
"""

Munge.spiking.nextandprev!(spikes)
perc = round(mean(spikes.neard .> .4)*100,sigdigits=2)


# =========PLOTS =======================================
Plot.setfolder( "isolated_ca1_rate_pfc")
histogram(log10.(spikes.neard),label="nearest spike")
vline!([log10(.400)], label="thresh", title="$perc % isolated")
Plot.save("nearest spike isolation")
# ======================================================


# ======================================================
#  _   _    _      _____           _ 
# | | | |  / \    |_   _| __ __ _ (_)
# | |_| | / _ \     | || '__/ _` || |
# |  _  |/ ___ \    | || | | (_| || |
# |_| |_/_/   \_\   |_||_|  \__,_|/ |
#                               |__/ 
# ======================================================
# ==============================
# GET VARIABLES OF INTEREST
# ==============================
per_hatraj = 
    get_isolation_summary(spikes, [:cuemem, :interneuron, :hatraj, :correct,
                                  :period])
dropmissing!(per_hatraj)
per_hatraj = 
    DIutils.arr.get_quantile_filtered(per_hatraj, :events_per_time, 0.02)
per_ctha =
    get_isolation_summary(spikes, [:area, :cuemem, :interneuron, :ha, :correct,
        :period])
dropmissing!(per_ctha)
# filter out top 0.1% of values
per_ctha = DIutils.arr.get_quantile_filtered(per_ctha, :events_per_time, 0.02)
# per_hatraj.cortsk = replace.(per_hatraj.cortsk, " " => "\n")

# ------------------------------
# plot ha + trajectory labels
# ------------------------------
iso_celltype_hatraj = @subset(sort(per_hatraj, 
                [:area, :interneuron, :cuemem, :hatraj]),
                first.(:hatraj) .!= '*')
cuemem_color = Dict(0=>:blue, 1=>:red)
correct_color = Dict(0=>:black, 1=>:green)
correct_fillstyle = Dict(0=>:/, 1=>nothing)
interneuron_style = Dict(true=>:dash, false=>:solid)
interneuron_string = Dict(true=>"Interneuron", false=>"Pyramidal")
ha_string = Dict('A'=>"Arena", 'H'=>"Home")
ca1_pfc_background_color = Dict("CA1"=>colorant"white", "PFC"=>colorant"lightgray")
ylabel_for_props = Dict(
    :events_per_time=>"Events per time",
    :isolated_mean=>"Fraction isolated spikes",
    :isolated_events_per_time=>"Isolated events per time"
)

"""
    subset_conservative_hatraj(per_hatraj)
Subset the data to only include the conservative hatraj labels:
CUE: A1,A2,H1,H2
MEM: A3,A4,H3,H4
"""
function subset_conservative_hatraj(per_hatraj)
    sub=@subset(per_hatraj, 
        (:cuemem .== 0 .&& in.(:hatraj, [("A1", "A2", "H1", "H2")])) .|
        (:cuemem .== 1 .&& in.(:hatraj, [("A3", "A4", "H3", "H4")]))
)
    @assert "H3" ∉ @subset(sub, :cuemem .== 0).hatraj .&&
            "H5" ∉ @subset(sub, :cuemem .== 1).hatraj .&&
            "A3" ∉ @subset(sub, :cuemem .== 0).hatraj .&&
            "A5" ∉ @subset(sub, :cuemem .== 1).hatraj
    sub
end

"""
    get_ylim_vars(per_ctha, prop)
Get the ylims for the plots.
"""
function get_ylim_vars(per_ctha, prop)
    area_interneuron_ylim = Dict((;area,interneuron) => 
            DIutils.plotutils.get_lims(@subset(per_ctha, :area .== area, :interneuron .== interneuron)[:,:events_per_time], 0.01)
        for (area, interneuron) in Iterators.product(["CA1", "PFC"], [true, false]))
    interneuron_area_ylim = Dict(
    (interneuron, area) => (0, maximum(@subset(per_ctha, :area .== area, :interneuron .==
            interneuron)[:,prop]))
        for (interneuron, area) in Iterators.product([true, false], 
                                                     ["CA1", "PFC"]))
    return area_interneuron_ylim, interneuron_area_ylim
end


# ------------------------------
# plot ha alone, no traj labels
# ------------------------------
"""
    plot_HAtraj(per_ctha, prop=:events_per_time)
Plot the number of events per time for each head angle.
# Arguments
- `per_ctha::DataFrame`: The data frame with the data to plot.
- `prop::Symbol`: The property to plot.
# Returns
# - `P::Array{Plots.Plot,1}`: The array of plots.
"""
function plot_HAtraj(iso_celltype_hatraj, prop=:events_per_time; 
    conservative_hatraj=true)
    if conservative_hatraj
        @info "Using conservative hatraj labels"
        iso_celltype_hatraj = subset_conservative_hatraj(iso_celltype_hatraj)
    end
    area_interneuron_ylim, interneuron_area_ylim = 
        get_ylim_vars(iso_celltype_hatraj, prop)
    P = []
    iters = Iterators.product([0, 1], ["CA1", "PFC"], [true, false])
    (cuemem, area, interneuron) = first(iters)
    for (cuemem, area, interneuron) in iters
        ich_per = sort(@subset(iso_celltype_hatraj, :interneuron .== interneuron, 
            :area .== area, :correct .== 1, :cuemem .== cuemem),
            [:cuemem, :hatraj])
        ich = combine(groupby(ich_per, [:cuemem, :hatraj]),
            prop=>x->nanmean(filter(x) do y
                !isnan(y) & !isinf(y) & !ismissing(y)
            end, x),
            prop=>(x->lower_stat_quant(mean, filter(x) do y
                    !isnan(y) & !isinf(y) & !ismissing(y)
            end, 0.005))=>:lower,
            prop=>(x->upper_stat_quant(mean, filter(x) do y
                !isnan(y) & !isinf(y) & !ismissing(y) 
            end, 0.975))=>:upper;
            renamecols=false)
        p= plot(ich.hatraj, ich[!,prop], xlabel="HATRAJ", 
            cue_color=cuemem_color[cuemem],
            ylabel="$(ylabel_for_props[prop])\n$(filt_desc[:all])", alpha=0.5,
            linewidth=2, label="",
            ribbon=(ich.lower, ich.upper), fillalpha=0.2, 
            fillcolor=cuemem_color[cuemem],
            title="$(interneuron_string[interneuron]),"*
                    "\nArea: $area, CueMem: $cuemem",
            linestyle=interneuron_style[interneuron], color=cuemem_color[cuemem],
            ylims=interneuron_area_ylim[(interneuron, area)])
        scatter!(ich.hatraj, ich[!,prop], 
            c=:black, label="", alpha=0.2, markersize=3)
        push!(P, p)
    end
    return P
end

P = plot_HAtraj(iso_celltype_hatraj, :events_per_time, conservative_hatraj=true)
plot(P..., layout=(4,2), size=(1000, 800))

P = plot_HAtraj(per_hatraj, :isolated_mean)
plot(P..., layout=(4,2), size=(1000, 800))

# ==============================
# HA labels only
# ==============================
"""
    plot_HA(per_ctha, prop=:events_per_time)
Plot the number of events per time for each head angle.
# Arguments
- `per_ctha::DataFrame`: The data frame with the data to plot.
- `prop::Symbol`: The property to plot.
# Returns
# - `P::Array{Plots.Plot,1}`: The array of plots.
"""
function plot_HA(per_ctha, prop)
    area_interneuron_ylim, interneuron_area_ylim = 
        get_ylim_vars(per_ctha, prop)
    ich = nothing
    iters = Iterators.product(['H','A'], ["CA1", "PFC"], [true, false])
    (ha, area, interneuron) = first(iters)
    P = OrderedDict()
    for (ha, area, interneuron) in iters
        key = (;ha, area, interneuron)
        ich_per = sort(@subset(per_ctha, 
            :interneuron .== interneuron, 
            :area .== area, :ha .== ha, :cuemem .!= -1, :correct .!= -1),
            [:cuemem, :ha])
        ich_per = ich_per[
        (!).(isnan.(ich_per[!,prop])) .&&
        (!).(isinf.(ich_per[!,prop])), :]
        ich = combine(groupby(ich_per, [:cortsk, :correct, :cuemem]), 
            prop=>x->nanmean(skipmissing(x)),
            prop=>(x->lower_stat_quant(mean,x,0.025))=>:lower,
            prop=>(x->upper_stat_quant(mean,x,0.975))=>:higher;
            # prop=>(x-> quantile(skipmissing(x), 0.025))=>:lower,
            # prop=>(x-> quantile(skipmissing(x), 0.975))=>:higher;
            renamecols=false)
        transform!(ich, :lower=>x->abs.(x), :higher=>x->abs.(x), renamecols=false)
            pval = nothing
            try
                T = HypothesisTests.UnequalVarianceTTest
                ttest = T(@subset(ich_per, :cuemem .== 0, :correct .== 1)[:, prop],
                          @subset(ich_per, :cuemem .== 1, :correct .== 1)[:, prop]
                )
                pval = round(pvalue(ttest), digits=3)
            catch
                pval = NaN
            end
            printstyled("pval: $pval\n", color=:red)
        p=plot()
        for (cuemem, correct) in Iterators.product([0, 1], [0, 1])
            XX = @subset(ich, :correct.== correct, :cuemem .== cuemem)
            if isempty(XX)
                continue
            end
            bar!(XX.cortsk, XX[!,prop], 
                yerror=(XX.lower, XX.higher),
                group=XX.correct,
                xlabel="CueMem", 
                alpha=0.5,
                title="$(interneuron_string[interneuron]), $area, $(ha_string[ha])",
                color=cuemem_color[cuemem],
                fillstyle=(println(correct_fillstyle[correct]); 
                                 correct_fillstyle[correct]),
                linewidth=1,
                fillcolor=cuemem_color[cuemem],
                label="",
            )
            scatter!(XX.cortsk, XX[!,prop], 
                yerror=(XX.lower, XX.higher),
                group=XX.correct,
                marker=:none, markeralpha=0.0,
                xlabel="CueMem", 
                alpha=0.5,
                title="$(interneuron_string[interneuron]), $area, $(ha_string[ha])",
                color=cuemem_color[cuemem],
                linewidth=5,
                msw=2,
                linestroke=(2,:dash),
                label="",
            )
            # println("ylabel_for_props[prop] = $(ylabel_for_props[prop])")
            Xper = @subset(ich_per, :correct.== correct, :cuemem .== cuemem)
            scatter!(Xper.cortsk, Xper[!,prop],
                group=Xper.correct,
                xlabel="CueMem", 
                ylabel="$(ylabel_for_props[prop])\n$(filt_desc[:all])", 
                alpha=0.1,
                title="$(interneuron_string[interneuron])," *
                    "\n $area, $(ha_string[ha])" *
                    "\nttest2(mem,cue)=\n$(round(pval, digits=3))",
                linestyle=interneuron_style[interneuron], 
                color=cuemem_color[cuemem],
                fillstyle=(println(correct_fillstyle[correct]); 
                                 correct_fillstyle[correct]),
                linewidth=0,
                fillcolor=cuemem_color[cuemem],
                markersize=2,
                label="",
                ylims=interneuron_area_ylim[(interneuron, area)],
            )
        end
        # Get scatter series elements
        sers = filter(p.series_list) do ser
            ser[:seriestype] == :scatter && length(ser[:y]) > 1
        end
        # Add jitter to the :x property of each sers element
        for ser in sers
            ser[:x] = ser[:x] .+ randn(length(ser[:x])) .* 0.1
        end
        p.attr[:background_color]         = ca1_pfc_background_color[area]
        p.attr[:foreground_color] = ca1_pfc_background_color[area]
        push!(P, key=>p)
    end
    return P
end

# ISSUE:
# 1. yaxis should be either quantile or K*mean, whichever comes first

Plot.setfolder("HA-split")
Plot.printstate()
#
# -------EVENTS PER TIME------------------------------------
    P = plot_HA(per_ctha, :events_per_time)
    pca1 = filter(P) do (k, v)
       k.area == "CA1"  
    end |> values |> collect
    plot(pca1..., layout=(1,4), size=(1000, 400), left_margin=10Plots.mm)
    Plot.save("HA_ca1.svg")
    ppfc = filter(P) do (k, v)
       k.area == "PFC"
    end |> values |> collect
    plot(ppfc..., layout=(1,4), size=(1000, 400), left_margin=10Plots.mm)
    Plot.save("HA_pfc.svg")
# ----------------------------------------------------------

# -----------ISOLATED MEAN----------------------------------
    P = plot_HA(per_ctha, :isolated_mean)
    pca1 = filter(P) do (k, v)
       k.area == "CA1"  
    end |> values |> collect
    plot(pca1..., layout=(1,4), size=(1000, 400), left_margin=10Plots.mm)
    Plot.save("HA_ca1_isolated.svg")
    ppfc = filter(P) do (k, v)
       k.area == "PFC"
    end |> values |> collect
    plot(ppfc..., layout=(1,4), size=(1000, 200), left_margin=10Plots.mm)
    Plot.save("HA_pfc_isolated.svg")
# ----------------------------------------------------------

# -----------ISOLATED MEAN----------------------------------
    P = plot_HA(per_ctha, :isolated_events_per_time)
    pca1 = filter(P) do (k, v)
       k.area == "CA1"  
    end |> values |> collect
    plot(pca1..., layout=(1,4), size=(1000, 400), left_margin=10Plots.mm)
    Plot.save("HA_ca1_isolated_events_per_time.svg")
    ppfc = filter(P) do (k, v)
       k.area == "PFC"
    end |> values |> collect
    plot(ppfc..., layout=(1,4), size=(1000, 400), left_margin=10Plots.mm)
    Plot.save("HA_pfc_isolated_events_per_time.svg")
# ----------------------------------------------------------
