
if !(:lfp in names(Main))
    include("./load_isolated.jl")
end

# Add a column to our spikes dataframe about its cell's meanrate
Load.register(cells, spikes, on="unit", transfer=["meanrate"])
# Add a behavioral info to spikes
Load.register(beh, spikes, on="time", transfer=["velVec", "period","correct"])
@assert :period ∈ propertynames(spikes)

# Acquire the isolation table
iso_sum = get_isolation_summary(spikes)
jldopen(filename, "a") do storage
    storage["iso_sum"] = iso_sum
end

# Obtain the isolated spikes dataframe
isolated = last(groupby(subset(spikes, :isolated=>x->(!).(isnan.(x))) ,
                            :isolated))
@assert all(isolated.isolated .== true)

# Select interneuron, roughly
spikes.interneuron = spikes.meanrate .> 5
histogram(cells.meanrate)
iso_sum_celltype = get_isolation_summary(spikes,[:cuemem, :interneuron])
sort!(iso_sum_celltype, [:area, :interneuron, :cuemem])
jldopen(filename, "a") do storage
    storage["iso_sum_celltype"] = iso_sum_celltype
end

iso_sum_celltype_per = get_isolation_summary(spikes, [:cuemem, :interneuron, :period, :correct])
@subset!(iso_sum_celltype_per, :events_per_time .!= Inf)
sort!(iso_sum_celltype_per, [:area, :interneuron, :cuemem, :period])
@subset!(iso_sum_celltype_per, (:cuemem .== -1 .&& :correct .== -1) .||
                               (:cuemem .== 0 .&& :correct .!= -1)  .||
                               (:cuemem .== 1 .&& :correct .!= -1))
subset!(iso_sum_celltype_per, :events_per_time => x-> (!).(isinf.(x)))
iso_sum_celltype_per
jldopen(filename, "a") do storage
    storage["iso_sum_celltype_per"] = iso_sum_celltype_per
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
                                          legendtitle="Task")
    Plot.save((;desc="time spent"))
    @df iso_sum bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:cuemem; kws..., 
                   legendtitle="Task")
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
    @df iso_sum_celltype bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:interneuron; kws..., legendtitle="is interneuron?", alpha=0.5)
    Plot.save((;desc="MUA per second, celltype"))
    @df iso_sum_celltype bar(:cuearea, :isolated_mean, ylabel="Isolated Spikes / All spikes\nsign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws..., legendtitle="is interneuron?", alpha=0.5)
    Plot.save((;desc="fraction of isolated spikes, celltype"))
    @df iso_sum_celltype bar(:cuearea, :isolated_events_per_time, ylabel="Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:interneuron; kws..., legendtitle="is interneuron?", alpha=0.5)
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
    sc=@df @subset(iso_sum_celltype_per, :interneuron .== false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :isolated_mean, group=:correct, alpha=0.5,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="Isolation",
                xlabel="task", ylabel="%percent isolated spikes",
                legend_title="correct/error", legend_position=:outerbottomright)
    end
    Plot.save((;desc="isolation", group=:correct, group_sep=true, ))

    @df iso_sum_celltype_per begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :isolated_mean, group=:correct, alpha=0.5)
    end

# PYRAMIDAL CELLS 🔺🔺
# begin

    # -------- isolated per second spikes ------------------
    @df @subset(iso_sum_celltype_per, :interneuron .== false) begin
        scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :isolated_events_per_time, group=:correct)
    end
    Plot.save((;desc="isolated-MUA-per-time", group=:correct, group_sep=true))
    @df @subset(iso_sum_celltype_per, :interneuron .== false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :isolated_events_per_time, group=:correct, alpha=0.5)
    end
    Plot.save((;desc="isolated-MUA-per-time", group=:interneuron, group_sep=true))

    # -------- MUA events per time -------------------------
    @df @subset(iso_sum_celltype_per, :interneuron .== false) begin
        scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), :events_per_time,
                group=:correct,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA
                events/time", xlabel="task", ylabel="events × time\$^{-1}\$,
                pyr cells", legend_title="correct/error", alpha=0.5,
                legend_position=:outerbottomright)
    end

    @df @subset(iso_sum_celltype_per, :interneuron .== false) begin
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
    @df @subset(iso_sum_celltype_per, :interneuron .== true) begin
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
    @df iso_sum_celltype_per begin
        scatter(:cuemem .+ (randn(size(:cuemem)) .* 0.1), 
                :events_per_time, group=:interneuron, alpha=0.5,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time",
                xlabel="task", ylabel="events × time\$^{-1}\$",
                legend_title="pyr/int", legend_position=:outerbottomright)
    end
    Plot.save((;desc="mua-per-time",group=:interneuron))

    @df iso_sum_celltype_per begin
        scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                :events_per_time, group=:interneuron, alpha=0.5,
                xtick=([-1,0,1],["nontask","cue","mem"]),title="MUA events/time", 
                xlabel="task", ylabel="events × time\$^{-1}\$",
                legend_title="Is Interneuron?",
                legend_position=:outerbottomright)
    end
    Plot.save((;desc="mua-per-time",group=:interneuron,group_sep=true))

    cpyr_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), 
                :isolated_mean, group=:correct, alpha=0.5, title="ca1 pyr",
                ylabel="isolated fraction", legendtitle="correct / error / nontask")
    end
    cint_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==true) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), 
                :isolated_mean, group=:correct, alpha=0.5, title="ca1 int",
                ylabel="isolated fraction", legendtitle="correct / error / nontask")
    end
    ppyr_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==false) begin
        scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), 
                :isolated_mean, group=:correct, alpha=0.5, title="pfc pyr",
                ylabel="isolated fraction", legendtitle="correct / error / nontask")
    end
    #pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
    pint_ce = Plot.create_blank_plot()

    plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, layout=grid(2,2))
    Plot.save((;desc="isolation cell type and correct"))

    p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0, :area .== "PFC") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
    end
    p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1, :area .== "PFC") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr mem")
    end
    p2 = p4 = Plot.create_blank_plot()
    plot(p1,p2,p3,p4, markersize=2)
    Plot.save((;desc="fract vs mua-per-sec, PFC"))

    p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
    end
    p2=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 0, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="int cue")
    end
    p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr mem")
    end
    p4=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 1, :area .== "CA1") begin
        scatter(:isolated_mean, :events_per_time, group=:correct, alpha=0.5, xlabel="fraction(isolated)",ylabel="multiunit/s", title="pyr cue")
    end
    plot(p1,p2,p3,p4, markersize=2)
    Plot.save((;desc="fract vs mua-per-sec, CA1"))
#end

# PYRAMIDAL 🔺🔺 and INTERNEURON  🔵
    # ⇛ summary of period-wise ISOLATION 
        @df iso_sum_celltype_per begin
            scatter(:cuemem .+ 0.25.*(:interneuron .- 0.5) .+ (randn(size(:cuemem)) .* 0.05),
                    :isolated_mean, group=:interneuron, alpha=0.5)
        end
        Plot.save((;desc="isolation", group=:interneuron, group_sep=true))
    # ➡ summary of period-wise EVENTS-PER-TIME 📈
    #begin
        cpyr_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==false) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time,
                    group=:correct, alpha=0.5, title="ca1 pyr", 
                    xlabel="task (nontask cue mem)", ylabel="events per time", legendtitle="correct error"
                   )
        end
        cint_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==true) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time, group=:correct, alpha=0.5, title="ca1 int",
                    xlabel="task (nontask cue mem)", ylabel="events per time", legendtitle="correct error")
        end
        ppyr_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==false) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :events_per_time, group=:correct, alpha=0.5, title="pfc pyr",
                    xlabel="task (nontask cue mem)", ylabel="events per time", legendtitle="correct error")
        end
        #pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
        pint_ce = Plot.create_blank_plot()
        plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, margin=30Plots.px, layout=grid(2,2), size=(1200,800))
    # end
    # ➡ summary of period-wise ISOLATION
        cint_ce = @df @subset(iso_sum_celltype_per, :area.=="CA1", :interneuron.==true) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :isolated_mean, group=:correct, alpha=0.5, title="ca1 int")
        end
        ppyr_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==false) begin
            boxplot(:cuemem .+ 0.25.*(:correct .- 0.5), :isolated_mean, group=:correct, alpha=0.5, title="pfc pyr")
        end
        #pint_ce = @df @subset(iso_sum_celltype_per, :area.=="PFC", :interneuron.==true) scatter(:cuemem .+ 0.25.*(:correct .- 0.5) .+ (randn(size(:cuemem)) .* 0.05), :isolated_mean, group=:correct, alpha=0.5, title="pfc int")
        pint_ce = Plot.create_blank_plot()
        plot(cpyr_ce, cint_ce, ppyr_ce, pint_ce, layout=grid(2,2))
        Plot.save((;desc="isolation cell type and correct"))

# Animations 🎥
    # -------- multi - metric : one predict other ? --------
    anim = @animate for i in 1:360
        p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                    group=:correct, alpha=0.5, camera=(i,30), 
                    title="pyr cue",
                    xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s")
        p2=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1) scatter(:isolated_mean, :isolated_events_per_time, :events_per_time,
                    group=:correct, alpha=0.5, camera=(i,30),
                    title="pyr mem",
                    xlabel="fraction(isolated)", ylabel="isolated multiunit/s", zlabel="multiunit/s")
        plot(p1,p2)
    end
    gif(anim, plotsdir(Plot.complete_folder_args..., "separatrix-three-vars-gif_pyr_facet=cuemem_group=correct.gif"); loop=1)

    anim = @animate for i in 1:360
        p1=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 1) begin
            scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr mem")
        end
        p2=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 1) begin
            scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30),ylabel="isolated multiunit/s", zlabel="multiunit/s", title="int mem")
        end
        p3=@df @subset(iso_sum_celltype_per, :interneuron .== false, :cuemem .== 0) begin
            scatter(:isolated_mean, :isolated_events_per_time, :events_per_time, group=:correct, alpha=0.5, camera=(i,30), xlabel="fraction(isolated)",ylabel="isolated multiunit/s", zlabel="multiunit/s", title="pyr cue")
        end
        p4=@df @subset(iso_sum_celltype_per, :interneuron .== true, :cuemem .== 0) begin
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


XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "CA1"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="isolated frac", title="CA1 pyr")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== true, :area .== "CA1"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="isolated frac", title="CA1 int")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "PFC"),
                     :cuemem), :isolated_mean => median, :isolated_mean => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="isolated frac", title="PFC pyr")
@df XX bar(:cuemem, :isolated_mean_median, yerror=:isolated_mean_function;
          kws...)
Plot.save(kws)

"""
summary of period mua
"""
XXg = groupby(@subset(iso_sum_celltype_per, :area .== "CA1"), :cuemem)
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

XXg = groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "CA1"), :cuemem)
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

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== true, :area .== "CA1"),
                     :cuemem), 
             :events_per_time => median, 
             :events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="MUA per sec", title="CA1 int MUA")
@df XX bar(:cuemem, :events_per_time_median, yerror=:events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "PFC"),
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

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "CA1"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
# TODO median, not mean, and bootstrap the median
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="CA1 pyr")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== true, :area .== "CA1"),
                     :cuemem), :isolated_events_per_time => median, :isolated_events_per_time => x->std(x)/sqrt(length(x))
            )
kws=(;xlabel="cuemem", ylabel="iso events per sec", title="CA1 int")
@df XX bar(:cuemem, :isolated_events_per_time_median, yerror=:isolated_events_per_time_function;
          kws...)
Plot.save(kws)

XX = combine(groupby(@subset(iso_sum_celltype_per, :interneuron .== false, :area .== "PFC"),
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

    svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
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

        svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
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

        svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
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

        svm_data = @subset(iso_sum_celltype_per, :correct .== 0 .|| :correct .== 1);
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


