quickactivate("/home/ryoung/Projects/goal-code/")
using Base.Threads: @spawn
using DataFrames
using StatsPlots
using Statistics
includet(srcdir("field.jl"))
includet(srcdir("filt.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("field/timeshift.jl"))
includet(srcdir("field/info.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("table.jl"))
spikes, beh = raw.load("RY16", 36, data_source=["spikes","behavior"])
props = ["x", "y"]
splitby=["unit", "area"]
kws=(;splitby, filters=merge(filt.speed_lib, filt.cellcount))
newkws = (; kws..., resolution=40, gaussian=2.3*0.5, props=props,
          filters=merge(kws.filters))

# MULTIPLE TIME SHIFT
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
@time place = timeshift.get_field_shift(beh, spikes, -1:0.01:1; 
                                        multithread=true, 
                                        postfunc=info.information, newkws...);
utils.pushover("Finished field shift")

# ---------
# 1 thread
# ---------
# (run 1) 108 seconds for 20 shifts without multi-threading, with 10% compile time
place_fineSecond = @spawn timeshift.get_field_shift(beh, spikes, -1:0.01:1;
                                        postfunc=info.information,
                                        multithread=false, newkws...);

place_broadTimes = @spawn @time timeshift.get_field_shift(beh, spikes, -4:0.05:4;
                                        postfunc=info.information,
                                        multithread=false, newkws...);
place_broadTimes = Dict(fetch(place_broadTimes)...)
utils.pushover("Finished broadtimes")


using StatsPlots
using DataFrames

function plot_shifts(place; desc="")
    descSave = replace(desc, ":"=>"", " "=>"-")
    # ALL CELLS
    df = table.to_dataframe(place, key_name="shift", name="info")
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
               title="$desc\nMI(Place)", xlabel="time", ylabel="cell"
    )
    vline!(p[1], [0], c=:white, linestyle=:dash, label="Zero lag")
    vline!(p[2], [0], c=:white, linestyle=:dash, legendposition=:none)
    savefig(plotsdir("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.pdf"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.png"))
    savefig(plotsdir("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.svg"))
end
plot_shifts(place_broadTimes, desc="Traj length sample: ")
