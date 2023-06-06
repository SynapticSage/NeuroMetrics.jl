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
        "filt"          => :all
    )

# @eval Main include(scriptsdir("isolated", "load_isolated.jl"))
@eval Main include("./load_isolated.jl")
@assert :isolated ∈ propertynames(spikes)
@assert beh !== nothing
Munge.nonlocal.setunfilteredbeh(beh)

using DI.Labels, DataFrames, DataFramesMeta, StatsPlots
using DIutils.plotutils
using GoalFetchAnalysis.Munge.isolated
using GoalFetchAnalysis.Munge.nonlocal

Plot.setparentfolder("isolated")
Plot.setappend("$animal-$day")

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
begin
    @df iso_ct bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:interneuron; kws..., legendtitle="is interneuron?", alpha=0.5)
    Plot.save((;desc="MUA per second, celltype"))
    @df iso_ct bar(:cuearea, :isolated_mean, ylabel="Isolated Spikes / All spikes\nsign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws..., legendtitle="is interneuron?", alpha=0.5)
    Plot.save((;desc="fraction of isolated spikes, celltype"))
    @df iso_ct bar(:cuearea, :isolated_events_per_time, ylabel="Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws..., legendtitle="is interneuron?", alpha=0.5)
    Plot.save((;desc="isolated spikes per second, celltype"))
end
# ======================================================

#    _  _     ____           _           _               _          
#  _| || |_  |  _ \ ___ _ __(_) ___   __| |    __      _(_)___  ___ 
# |_  ..  _| | |_) / _ \ '__| |/ _ \ / _` |____\ \ /\ / / / __|/ _ \
# |_      _| |  __/  __/ |  | | (_) | (_| |_____\ V  V /| \__ \  __/
#   |_||_|   |_|   \___|_|  |_|\___/ \__,_|      \_/\_/ |_|___/\___|
#                                                                   

# =========PLOTS =======================================
# -------- mua per second ------------------------------
Plot.setfolder("MUA and isolated MUA", "period-wise")

# CELL AGNOSTIC ⬜
    # -------- percent isolated spikes ---------------------
    sc=@df @subset(per_ct, :interneuron .== false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :isolated_mean, group=:correct, alpha=0.5,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="Isolation",
                xlabel="task", ylabel="%percent isolated spikes",
                legend_title="correct/error", legend_position=:outerbottomright)
    end
    Plot.save((;desc="isolation", group=:correct, group_sep=true, ))

    @df per_ct begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :isolated_mean, group=:correct, alpha=0.5)
    end

# PYRAMIDAL CELLS 🔺🔺
# begin

    # -------- isolated per second spikes ------------------
    @df @subset(per_ct, :interneuron .== false) begin
        scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :isolated_events_per_time, group=:correct)
    end
    Plot.save((;desc="isolated-MUA-per-time", group=:correct, group_sep=true))
    @df @subset(per_ct, :interneuron .== false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :isolated_events_per_time, group=:correct, alpha=0.5)
    end
    Plot.save((;desc="isolated-MUA-per-time", group=:interneuron, group_sep=true))

    # -------- MUA events per time -------------------------
    @df @subset(per_ct, :interneuron .== false) begin
        scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :events_per_time,
                group=:correct,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA
                events/time", xlabel="task", ylabel="events × time\$^{-1}\$,
                pyr cells", legend_title="correct/error", alpha=0.5,
                legend_position=:outerbottomright)
    end

    @df @subset(per_ct, :interneuron .== false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), 
                :events_per_time, group=:correct, alpha=0.5,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA
                events/time", xlabel="task", ylabel="events × time\$^{-1}\$,
                pyr cells", legend_title="correct/error",
                legend_position=:outerbottomright)
    end
    Plot.save((;desc="mua-per-time",group=:correct,pyr=true))

# end

# INTERNEURONS 🔵
#begin
    @df @subset(per_ct, :interneuron .== true) begin
        scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1),
                :events_per_time, group=:correct,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA
                events/time", xlabel="task", ylabel="events × time\$^{-1}\$,
                pyr cells", legend_title="correct/error",
                legend_position=:outerbottomright)
    end
    Plot.save((;desc="mua-per-time",group=:correct,pyr=true,sep_groups=true))
#end

# PYRAMIDAL 🔺🔺 and INTERNEURON  🔵
# ➡ each period sample
#begin
    @df per_ct begin
        scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), 
                :events_per_time, group=:interneuron, alpha=0.5,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time",
                xlabel="task", ylabel="events × time\$^{-1}\$",
                legend_title="pyr/int", legend_position=:outerbottomright)
    end
    Plot.save((;desc="mua-per-time",group=:interneuron))

    @df per_ct begin
        scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :events_per_time, group=:interneuron, alpha=0.5,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", 
                xlabel="task", ylabel="events × time\$^{-1}\$",
                legend_title="Is Interneuron?",
                legend_position=:outerbottomright)
    end
    Plot.save((;desc="mua-per-time",group=:interneuron,group_sep=true))

    cpyr_ce = @df @subset(per_ct, :area.=="CA1", :interneuron.==false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), 
                :isolated_mean, group=:correct, alpha=0.5, title="ca1 pyr",
                ylabel="isolated fraction", legendtitle="correct / error / nontask")
    end
    cint_ce = @df @subset(per_ct, :area.=="CA1", :interneuron.==true) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), 
                :isolated_mean, group=:correct, alpha=0.5, title="ca1 int",
                ylabel="isolated fraction", legendtitle="correct / error / nontask")
    end
    ppyr_ce = @df @subset(per_ct, :area.=="PFC", :interneuron.==false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), 
                :isolated_mean, group=:correct, alpha=0.5, title="pfc pyr",
                ylabel="isolated fraction", legendtitle="correct / error / nontask")
    end
    #pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
    pint_ce = Plot.create_blank_plot()

    plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, layout=grid(2,2))
    Plot.save((;desc="isolation cell type and correct"))

    p1=@df @subset(per_ct, :interneuron .== false, :cuemem .== 0, :area .== "PFC") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
    end
    p3=@df @subset(per_ct, :interneuron .== false, :cuemem .== 1, :area .== "PFC") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr mem")
    end
    p2 = p4 = Plot.create_blank_plot()
    plot(p1,p2,p3,p4, markersize=2)
    Plot.save((;desc="fract vs mua-per-sec, PFC"))

    p1=@df @subset(per_ct, :interneuron .== false, :cuemem .== 0, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
    end
    p2=@df @subset(per_ct, :interneuron .== true, :cuemem .== 0, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="int cue")
    end
    p3=@df @subset(per_ct, :interneuron .== false, :cuemem .== 1, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr mem")
    end
    p4=@df @subset(per_ct, :interneuron .== true, :cuemem .== 1, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
    end
    plot(p1,p2,p3,p4, markersize=2)
    Plot.save((;desc="fract vs mua-per-sec, CA1"))
#end

# PYRAMIDAL 🔺🔺 and INTERNEURON  🔵
    # ⇛ summary of period-wise ISOLATION 
        @df per_ct begin
            scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                    :isolated_mean, group=:interneuron, alpha=0.5)
        end
        Plot.save((;desc="isolation", group=:interneuron, group_sep=true))
    # ➡ summary of period-wise EVENTS-PER-TIME 📈
    #begin
        cpyr_ce = @df @subset(per_ct, :area.=="CA1", :interneuron.==false) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time,
                    group=:correct, alpha=0.5, title="ca1 pyr", 
                    xlabel="task (nontask cue mem)", ylabel="events per time", legendtitle="correct error"
                   )
        end
        cint_ce = @df @subset(per_ct, :area.=="CA1", :interneuron.==true) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time, group=:correct, alpha=0.5, title="ca1 int",
                    xlabel="task (nontask cue mem)", ylabel="events per time", legendtitle="correct error")
        end
        ppyr_ce = @df @subset(per_ct, :area.=="PFC", :interneuron.==false) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time, group=:correct, alpha=0.5, title="pfc pyr",
                    xlabel="task (nontask cue mem)", ylabel="events per time", legendtitle="correct error")
        end
        #pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
        pint_ce = Plot.create_blank_plot()
        plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, margin=30Plots.px, layout=grid(2,2), size=(1200,800))
    # end
    # ➡ summary of period-wise ISOLATION
        cint_ce = @df @subset(per_ct, :area.=="CA1", :interneuron.==true) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :isolated_mean, group=:correct, alpha=0.5, title="ca1 int")
        end
        ppyr_ce = @df @subset(per_ct, :area.=="PFC", :interneuron.==false) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :isolated_mean, group=:correct, alpha=0.5, title="pfc pyr")
        end
        #pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
        pint_ce = Plot.create_blank_plot()
        plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, layout=grid(2,2))
        Plot.save((;desc="isolation cell type and correct"))

# Animations 🎥
    # -------- multi - metric : one predict other ? --------
    anim = @animate for i in 1:360
        p1=@df @subset(per_ct, :interneuron .== false, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                    group=:correct, alpha=0.5, camera=(i,30), 
                    title="pyr cue",
                    xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s")
        p2=@df @subset(per_ct, :interneuron .== false, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                    group=:correct, alpha=0.5, camera=(i,30),
                    title="pyr mem",
                    xlabel="fraction(isolated)", ylabel="isolated multiunit/s", zlabel="multiunit/s")
        plot(p1,p2)
    end
    gif(anim, plotsdir(Plot.complete_folder_args..., "separatrix-three-vars-gif_pyr_facet=cuemem_group=correct.gif"); loop=1)

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
        plot(p3,p4,p1,p2; layout=grid(2,2))
    end
    gif(anim, plotsdir(Plot.complete_folder_args..., "separatrix-three-vars-gif_facet=pyrint,cuemem_group=correct.gif"); loop=1)

# ======================================================

"""
===========================
summary of period iso_mean
===========================
"""
Plot.setfolder("MUA and isolated MUA", "period-wise-summary")


XX = combine(groupby(@subset(per_ct, :interneuron .== false, :area .== "CA1"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="isolated frac", title="CA1 pyr")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(per_ct, :interneuron .== true, :area .== "CA1"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="isolated frac", title="CA1 int")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(per_ct, :interneuron .== false, :area .== "PFC"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="isolated frac", title="PFC pyr")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

"""
summary of period mua
"""
XXg = groupby(@subset(per_ct, :area .== "CA1"), :cuemem)
XX = combine(XXg, :events_per_time => median, 
             :events_per_time => x->nanstd(x)./sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="CA1 MUA per second")
@df XX bar(:cuemem, :events_per_time_median, yerr=:events_per_time_function;
          linewidth=2, 
          kws...)
@df XX plot!(:cuemem, :events_per_time_median, yerr=:events_per_time_function, linestyle=Symbol("none"))
Plot.save(kws)

XXg = groupby(@subset(per_ct, :interneuron .== false, :area .== "CA1"), :cuemem)
XX = combine(XXg, :events_per_time => median, 
             :events_per_time => x->nanstd(x)./sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="CA1 pyr MUA")
@df XX bar(:cuemem, :events_per_time_median, yerr=:events_per_time_function;
          linewidth=2, 
          kws...)
@df XX plot!(:cuemem, :events_per_time_median, yerr=:events_per_time_function, linestyle=:none)
Plot.save(kws)

XX = combine(groupby(@subset(per_ct, :interneuron .== true, :area .== "CA1"),
                     :cuemem), 
             :events_per_time => median, 
             :events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="CA1 int MUA")
@df XX bar(:cuemem, :events_per_time_median, yerror=:events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(per_ct, :interneuron .== false, :area .== "PFC"),
                     :cuemem), :events_per_time => median, :events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="PFC pyr")
@df XX bar(:cuemem, :events_per_time_median, yerror=:events_per_time_function;
          kws...)
Plot.save(kws)

"""
===========================
summary of iso events per time
===========================
"""

XX = combine(groupby(@subset(per_ct, :interneuron .== false, :area .== "CA1"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="CA1 pyr")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(per_ct, :interneuron .== true, :area .== "CA1"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="CA1 int")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(per_ct, :interneuron .== false, :area .== "PFC"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="PFC pyr")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)


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
