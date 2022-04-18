quickactivate("/home/ryoung/Projects/goal-code/")
using Base.Threads: @spawn
using DataFrames
using StatsPlots
using Statistics
using StatsPlots
using DataFrames
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("field/timeshift.jl"))
includet(srcdir("field/info.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("table.jl"))
spikes, beh = raw.load("RY16", 36, data_source=["spikes","behavior"])

# MULTIPLE TIME SHIFT NOTES
#
# -----------
# 16 threads
# -----------
# (run 1) 400 seconds for 20, with multi-threading, with 90% of that compile time
# (run 2) 180 seconds 85% compile time
# (run 3) 131 seconds 76% compile time
# (run 4) 131 seconds 76% compile time
#
# ---------
# 4 threads
# ---------
#
# ---------
# 1 thread
# ---------
# (run 1) 108 seconds for 20 shifts without multi-threading, with 10% compile time
#
# ------------
# More shifts
# ------------
# For 200 shifts, we're in the domain of 30 minutes
#
#
# A (( ---- PLACE ---- ))
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount))
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters))
place = @spawn @time timeshift.get_field_shift(beh, spikes, -2:0.1:2; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
place_fineNarrow = @spawn @time timeshift.get_field_shift(beh, spikes, -1:0.01:1;
                                        postfunc=info.information,
                                        multithread=false, newkws...);
place_broadTimes = @spawn @time timeshift.get_field_shift(beh, spikes, -4:0.05:4;
                                        postfunc=info.information,
                                        multithread=false, newkws...);
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount))
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters, filt.correct))
place_correctOnly = @time timeshift.get_field_shift(beh, spikes, -2:0.1:2;
                                                    multithread=false,
                                                    postfunc=info.information, newkws...);
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters, filt.incorrect))
place_incorrectOnly = @spawn @time timeshift.get_field_shift(beh, spikes, -2:0.1:2;
                                                    multithread=false,
                                                    postfunc=info.information, newkws...);
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters, filt.nontask))
place_nontaskOnly = @spawn @time timeshift.get_field_shift(beh, spikes, -2:0.1:2;
                                                    multithread=false,
                                                    postfunc=info.information, newkws...);

newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters, filt.cue, filt.correct))
place_correctCue = @spawn @time timeshift.get_field_shift(beh, spikes, -2:0.1:2;
                                                    multithread=false,
                                                    postfunc=info.information, newkws...);
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters, filt.mem, filt.correct))
place_correctMem = @spawn @time timeshift.get_field_shift(beh, spikes, -2:0.1:2;
                                                    multithread=false,
                                                    postfunc=info.information, newkws...);


# B (( ---- GOAL ---- ))
filters = merge(kws.filters,
                filt.correct,
                filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150))
props = ["currentAngle", "currentPathLength"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=40, gaussian=2.3*0.5, props=props)
goal = @spawn @time timeshift.get_field_shift(beh, spikes, -1:0.05:1; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
goal_fineNarrow = @spawn @time timeshift.get_field_shift(beh, spikes, -1:0.01:1; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
goal_broadTimes = @spawn @time timeshift.get_field_shift(beh, spikes, -4:0.05:4; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
# C (( ---- SPECGOAL ---- ))
props = ["currentAngle", "currentPathLength", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 5], gaussian=2.3*0.5, props=props)
specgoal = @spawn @time timeshift.get_field_shift(beh, spikes, -1:0.05:1; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
# C (( ---- SPECPLACE ---- ))
props = ["x", "y", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 5], gaussian=2.3*0.5, props=props)
specplace = @spawn @time timeshift.get_field_shift(beh, spikes, -1:0.05:1; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
# D (( ---- FULL ---- ))
props = ["x", "y", "currentAngle", "currentPathLength", "stopWell"]
newkws = (; kws..., filters=merge(kws.filters, filters),
          resolution=[40, 40, 40, 40, 5], gaussian=2.3*0.5, props=props)
full = @spawn @time timeshift.get_field_shift(beh, spikes, -1:0.05:1; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);


# (( PLACE ))
if isdefined(Main, :place_broadTimes)
    place_broadTimes = Dict(fetch(place_broadTimes)...)
end
if isdefined(Main, :place_fineNarrow)
    place_fineNarrow = Dict(fetch(place_fineNarrow)...)
end
if isdefined(Main, :place)
    place = Dict(fetch(place)...)
end
# (( GOAL ))
if isdefined(Main, :goal)
    goal = Dict(fetch(goal)...)
end
if isdefined(Main, :goal_fineNarrow)
    goal_fineNarrow = Dict(fetch(goal_fineNarrow)...)
end
if isdefined(Main, :goal_broadTimes)
    goal_broadTimes = Dict(fetch(goal_broadTimes)...)
end
# (( SPECGOAL ))
if isdefined(Main, :specgoal)
    specgoal = Dict(fetch(specgoal)...)
end
# (( SPECPLACE ))
if isdefined(Main, :specplace)
    specplace = Dict(fetch(specplace)...)
end
if isdefined(Main, :place_correctOnly)
    place_correctOnly = Dict(fetch(place_correctOnly)...)
end
# (( FULL ))
if isdefined(Main, :full)
    full = Dict(fetch(full)...)
end
utils.pushover("Finished spawned shift processes")


function plot_shifts(place; desc="")

    descSave = replace(desc, ":"=>"", " "=>"-")
    

    # DENSITY ALL CELLS
    df = table.to_dataframe(place, key_name="shift", name="info")
    df.shift = df.shift .* -1; # beh.time-shift... negative shift is future, so correcting this
    df = sort(df, [:area,:unit,:shift])
    taus = unique(sort(df.shift))
    df_imax = combine(groupby(df, [:unit, :area]), 
                      :info=>argmax, 
                      :info=>(x->taus[argmax(x)])=>:bestTau)
    df_imax = df_imax[df_imax.bestTau.!=taus[1],:] # Excluding samples with the first tau, because that's the null condition for no variation

    @df df_imax density(:bestTau, group=:area, 
                 title="$desc BestTau(Information)", xlabel="Seconds", ylabel="Density")
    savefig(plotsdir("fields", "shifts", "$(descSave)_density_x=seconds,y=densBestTau_by=area.svg"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_density_x=seconds,y=densBestTau_by=area.png"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_density_x=seconds,y=densBestTau_by=area.pdf"))
    @df df_imax histogram(:bestTau, group=:area, 
                 title="$desc BestTau(Information)", xlabel="Seconds", ylabel="Density")
    savefig(plotsdir("fields", "shifts", "$(descSave)_histogram_x=seconds,y=densBestTau_by=area.svg"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_histogram_x=seconds,y=densBestTau_by=area.png"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_histogram_x=seconds,y=densBestTau_by=area.pdf"))

    # HISTOGRAM ALL CELLS
    df_m = combine(groupby(df, [:area, :shift]), :info=>mean)
    @df df_m bar(:shift, :info_mean, group=:area, 
                 title="$desc Median(Information)", xlabel="Seconds", ylabel="Shannon MI Bits")
    savefig(plotsdir("fields", "shifts", "$(descSave)_histogram_x=shift,y=medianinfo_by=area.svg"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_histogram_x=shift,y=medianinfo_by=area.png"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_histogram_x=shift,y=medianinfo_by=area.pdf"))

    # PER CELL
    df = sort(df, [:shift,:area,:unit])
    df_u = sort(unstack(df, :shift, :info), [:area, :unit])
    shifts = parse.(Float32,names(df_u)[3:end])
    units = df_u.unit
    areas = df_u.area
    area_divide = findfirst([diff(areas .== "PFC"); 0].==1)
    exclude = Not([x for x in [:area,:unit,:taus] if String(x) in names(df_u)])
    df_u.taus = vec([x[2] for x in argmax(Matrix(df_u[!,exclude]),dims=2)])
    df_u = sort(df_u, [:area, :taus])
    function get_area(area) 
        M = Matrix(df_u[df_u.area.==area, Not([:area,:unit, :taus])])
        M = replace(M, missing=>NaN)
        M[findall(vec([x[1] for x in any(M.!=0,dims=2)])),:]
    end
    p = Plots.plot(
               heatmap(shifts, 1:size(get_area("CA1"),1), get_area("CA1"), 
                clims=(0,20), colorbar_title="MI (bits)", colorbar_titlefontrotation=0),
               heatmap(shifts, 1:size(get_area("PFC"),1), get_area("PFC"), 
                clims=(0,8), colorbar_title="MI (bits)", colorbar_titlefontrotation=0),
               title="$desc\nMI", xlabel="time", ylabel="cell"
    )
    vline!(p[1], [0], c=:white, linestyle=:dash, label="Zero lag")
    vline!(p[2], [0], c=:white, linestyle=:dash, legendposition=:none)
    savefig(plotsdir("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.pdf"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.png"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.svg"))
end

plot_shifts(place_broadTimes,    desc="Traj length sample: ")
plot_shifts(place_fineNarrow,    desc="PLACE_HIGHRES")
plot_shifts(place,               desc="PLACE_LOWRES")
plot_shifts(goal_fineNarrow,     desc="GOAL_HIGHRES")
plot_shifts(goal_broadTimes,     desc="GOAL_TRAJ_LEN")
plot_shifts(specgoal,            desc="SPECGOAL")
plot_shifts(full,                desc="FULL")
plot_shifts(specplace,           desc="SPECPLACE")
plot_shifts(place_correctOnly,   desc="PLACE_CORRECTONLY")
plot_shifts(place_incorrectOnly, desc="PLACE_ERRORONLY")
plot_shifts(place_nontaskOnly,   desc="PLACE_NONTASK")
plot_shifts(place_correctCue,    desc="PLACE_CUEcorrect")
plot_shifts(place_correctMem,    desc="PLACE_MEMcorrect")
