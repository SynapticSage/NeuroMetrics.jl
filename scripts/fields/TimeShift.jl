quickactivate("/home/ryoung/Projects/goal-code/")
@time include(scriptsdir("fields", "Include.jl"))
@time include(scriptsdir("fields", "Initialize.jl"))
includet(srcdir("shuffle.jl"))
includet(srcdir("field/info.jl"))
includet(srcdir("field/timeshift.jl"))
include(scriptsdir("fields", "TimeShift_setsOfInterest.jl"))
_, spikes = raw.register(beh, spikes; transfer=["velVec"], on="time")
import Base.Threads: @spawn
using ThreadSafeDicts
using .timeshift
using Combinatorics: powerset
sp = copy(spikes)
convertToMinutes = false

# ----------
# ----------
# ----------
# Parameters
# ----------
PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)


#-------------------------------------------------------
"""
Translate into shorcut names
"""
𝕄(items)  = [replace(item, var_shortcut_names...) for item in items]
"""
UnTranslate from shorcut names
"""
𝕄̅(items)  = [replace(item, Dict(kv[2]=>kv[1] for kv in shortcut_names)...)
             for item in items]
sz(items) = [IDEALSIZE[item] for item in items]
#-------------------------------------------------------

splitby=["unit", "area"]
filters = merge(filt.notnan("currentAngle"), 
                filt.minmax("currentPathLength", 2, 150),
                filt.speed_lib, filt.cellcount)
gaussian=2.3*0.5

# COMPUTE INFORMATION @ DELAYS
I = Dict()
S = Dict()
for props ∈ marginals_highprior
    marginal = 𝕄(props)
    newkws = (; kws..., filters, splitby, gaussian, props,
              resolution=sz(props), multi=:single, postfunc=info.information)
    I[marginal] = @spawn @time timeshift.get_field_shift(beh, spikes, -2:0.1:2; 
                                                         newkws...)
    S[marginal] = @spawn @time timeshift.get_field_shift_shufflesfield_shift(beh, spikes, -2:0.1:2; 
                                                         newkws...)
end
for (key, value) in I
    I[key] = fetch(I[key])
end

# GET FIELDS AT BEST-(𝛕)
F = Dict()
for (key, value) in I
    newkws = (;kws..., filters, splitby, gaussian, props, resolution=sz(props))
    F      = timeshift.fetch_best_fields(I[key], beh, spikes; newkws...)
end


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
plot_shifts(place_broadTimes, desc="Traj length sample: ")
plot_shifts(place_fineNarrow, desc="PLACE_HIGHRES")
plot_shifts(place,           desc="PLACE_LOWRES")
plot_shifts(goal_fineNarrow, desc="GOAL_HIGHRES")
plot_shifts(goal_broadTimes, desc="GOAL_TRAJ_LEN")
plot_shifts(specgoal, desc="SPECGOAL")
plot_shifts(full,     desc="FULL")
plot_shifts(specplace,desc="SPECPLACE")
