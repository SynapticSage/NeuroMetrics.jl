
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using DataFrames
using Plots, Measures
using KernelDensity, Distributions
using Gadfly, Cairo, Fontconfig
import ColorSchemes
includet(srcdir("raw.jl"))
includet(srcdir("field.jl"))

beh, ripples = raw.load("RY16", 36, data_source=["behavior","ripples"]);
@assert ("x" ∈ names(beh)) "Fuck"
props = ["cuemem", "stopWell", "startWell", "x", "y"]
speedfilters = Dict("velVec"=>
               x->abs.(x) .≤ 5)
add_col = [[key for (key,value) in speedfilters]..., props...]
lookupcols = Dict((source=1,target=2) => add_col)
still_beh, ripples = raw.filterTables(beh, ripples; filters=speedfilters,
                                lookupcols=lookupcols)
ripples.duration = ripples.stop.-ripples.start;
ripples.dur= ripples.stop.-ripples.start;
ripples.ampdur = ripples.duration .* ripples.amp;
pfcripples = ripples[ripples.area.=="PFC",:]
ripples = ripples[ripples.area.=="CA1",:]

L(x) = Dict(-1=>"NonTask",
             0=> "Cue",
             1=>"Memory")[x]

function gsave(p, plotpath)
    p |> PNG(plotpath * ".png"); 
    p |> SVG(plotpath * ".svg");
    return nothing
end


# Amplitude and occ norm ripples and pfc ripples
filters = Dict("velVec"=>x->abs.(x) .≤ 5)
x = field.get_fields(beh, ripples, filters=speedfilters, resolution=100,
                     gaussian=2.3, behfilter=false);

filters = Dict("velVec"=>x->abs.(x) .≤ 5)
y = field.get_fields(beh, pfcripples, filters=speedfilters, resolution=100, gaussian=2.3,
                                          behfilter=false);

hca1 = heatmap(x.hist, colorbar=false, label="CA1", title="A",
               titlelocation=:left, clim=(0,0.008), framestyle=:none)
hpfc = heatmap(y.hist, colorbar=false, title="D", titlelocation=:right, clim=(0,0.005), framestyle=:none)

ampca1 = Plots.histogram(ripples.amp, label="CA1 global ripple dist",
                         title="B", titlelocation=:left, xlabel="Amplitude")
vline!([mean(ripples.amp)], linewidth=2, c="red", label="mean")
vline!([median(ripples.amp)], linewidth=2, c="black", linestyle=:dot,
       label="median")
amppfc = Plots.histogram(pfcripples.amp, c="skyblue", label="PFC global ripple
                         dist", title="E", titlelocation=:right, xlabel="Amplitude")
vline!([mean(pfcripples.amp)], linewidth=2, c="red", label="mean")
vline!([median(pfcripples.amp)], linewidth=2, c="black", linestyle=:dot,
       label="median")
xlims!(0, maximum(ripples.amp))
ylims!(0, 800)

durca1 = Plots.histogram(ripples.dur, label="CA1 global ripple dist",
                         title="C", titlelocation=:left, xlabel="Duration")
vline!([mean(ripples.dur)], linewidth=2, c="red", label="mean")
vline!([median(ripples.dur)], linewidth=2, c="black", linestyle=:dot,
       label="median")
xlims!(0, 1)
ylims!(0, 800)
durpfc = Plots.histogram(pfcripples.dur, c="skyblue", label="pfc global ripple
                         dist", title="F", titlelocation=:right, xlabel="duration")
vline!([mean(pfcripples.dur)], linewidth=2, c="red", label="mean")
vline!([median(pfcripples.dur)], linewidth=2, c="black", linestyle=:dot,
       label="median")
xlims!(0, 1)
ylims!(0, 800)

h=Plots.plot(hca1, hpfc, ampca1, amppfc, durca1, durpfc, layout=grid(3,2))
h.attr[:dpi]=200
h.attr[:size] = (600,1000)
savefig(h, plotsdir("ripples", "panelDistributions.pdf"))
savefig(h, plotsdir("ripples", "panelDistributions.svg"))

# Splitby cuemem

filters = Dict("velVec"=>x->abs.(x) .≤ 5)
x = field.get_fields(beh, ripples, splitby=["cuemem", "correct"], filters=speedfilters, 
                                          behfilter=false);

filters = Dict("velVec"=>x->abs.(x) .≤ 5)
y = field.get_fields(beh, pfcripples, splitby=["cuemem"], filters=speedfilters, 
                                          behfilter=false);


# Amplitudes
Gadfly.push_theme(:dark)
set_default_plot_size(15cm,10cm)
maze_coord = Coord.cartesian(xmin=60, ymin=12, ymax=85)

# Distribution
behaviorLayer = Gadfly.plot(beh, x=:x, y=:y, alpha=[0.1], 
                            maze_coord, Geom.hexbin(xbincount=120,
                                                    ybincount=120),
                            Guide.title("Behavior Distribution :: V ≤ 0.5")
                           );
rippleLayer = Gadfly.plot(ripples, maze_coord, x=:x, y=:y, 
            Geom.hexbin(xbincount=120, ybincount=120),
            Guide.title("Ripple Distribution")
           );
set_default_plot_size(40cm,15cm)
distributionRippleAndBehavior = hstack(behaviorLayer, rippleLayer)
plotpath = plotsdir("ripples","ripple_and_stillness_distribution")
gsave(distributionRippleAndBehavior, plotpath)


# By cuemem
coord = Coord.cartesian(ymax=50)
scale = Scale.color_discrete_hue(levels=[-1,0,1])
cuemem_and_error_amplitude = Gadfly.plot(ripples, coord,
                                         x=:amp, color=:cuemem, alpha=[0.5], scale,
                                         Geom.histogram(position=:identity, bincount=50));
cuemem_amplitude = Gadfly.plot(ripples[findall(ripples.cuemem.>=0),:], coord,
                               x=:amp, color=:cuemem, alpha=[0.5], scale,
            Geom.histogram(bincount=50, position=:identity));
ripple_amps_cuemem = hstack(cuemem_and_error_amplitude, cuemem_amplitude);
plotpath = plotsdir("ripples","ripple_amps_cuemem")
gsave(distributionRippleAndBehavior, plotpath)


# By cuemem
coord = Coord.cartesian(ymax=50)
scale = Scale.color_discrete_hue(levels=[-1,0,1])
cuemem_and_error_amplitude = Gadfly.plot(ripples, coord,
                                         x=:amp, color=:cuemem, alpha=[0.5], scale,
                                         Geom.histogram(position=:identity, bincount=50));
cuemem_amplitude = Gadfly.plot(ripples[findall(ripples.cuemem.>=0),:], coord,
                               x=:amp, color=:cuemem, alpha=[0.5], scale,
            Geom.histogram(bincount=50, position=:identity));
ripple_amps_cuemem = hstack(cuemem_and_error_amplitude, cuemem_amplitude);


# Barsum cuemem
function boot(dat, func, n_boot=200)
    bs = Tuple(rand(1:size(dat,1), size(dat,1)) for i in 1:n_boot)
    D = [func(dat[bs[i],:]) for i in 1:n_boot]
    return D
end
function mean_and_ci(D, col, alpha=0.99)
    col = String(col)
    dresult = D[1]
    res = [ d[!, col] for d in D] 
    res = hcat(res...)
    mres = mean(res;dims=2)
    lower = vcat([quantile(r, alpha/2) for r in eachslice(res,dims=1)]...)
    upper = vcat([quantile(r, 1 - alpha/2) for r in eachslice(res,dims=1)]...)
    dresult[!, col] = convert.(Float64, dresult[!,col])
    dresult[!,col] = vec(mres)
    dresult[!,"l"*col] = lower;
    dresult[!,"h"*col] = upper;
    return dresult
end
function measure(R, groups="cuemem"; stat=:count, what=nothing, beh=beh, occNorm=true)
    if what ≠ nothing
        yval  = Symbol(String(what)*"_"*String(stat))
        mapping = String(what) => getfield(Main, stat) => yval
    else
        mapping = nrow => :count
    end
    func(dat) = combine(groupby(dat, groups, sort=true), mapping)
    if stat == :count
        yval  = :count
        X = boot(R, func)
        X = mean_and_ci(X, yval)
    else
        X = boot(R, func)
        X = mean_and_ci(X, yval)
    end
    normalize_cols = [String(yval), "l"*String(yval), "h"*String(yval)]
    if occNorm
        X = table.occupancy_normalize(X, beh, groups, normalize_cols=normalize_cols)
    end
    X = transform(X, :cuemem => (x->L.(x)) => :cm)
    return X, (yval, "l"*String(yval), "h"*String(yval))
end
function even_out(args...)
    m = maximum([length(arg) for arg in args])
    for i in 1:length(args)
        while length(args[i]) < m
            push!(args[i], Gadfly.plot())
        end
    end
    return args
end

scalecm=Scale.color_discrete_manual(ColorSchemes.PRGn_3[[1,3]]...;
                                    levels=[0,1])
filt(dat) = dat[dat.cuemem .≥ 0 .&& dat.stopWell.>0, :]
R = filt(ripples)
B = filt(beh)
V1, V2, V3 = [], [], []
X, (yval, ymin, ymax) = measure(R, ["cuemem"], beh=B);
x =Gadfly.plot(X, scalecm, color=:cuemem, x=:cm, y=yval, ymin=ymin, ymax=ymax,
               Guide.ylabel("Ripple count\n(occnorm)", orientation=:vertical),
               Geom.bar, Geom.errorbar);
push!(V1, x);
# Barsum cuemem x goals
X, (yval, ymin, ymax) = measure(R, ["stopWell", "cuemem"], beh=B);
x = Gadfly.plot(X, x=:cm, color=:stopWell, y=yval, ymin=ymin, ymax=ymax,
                Geom.bar, Geom.errorbar);
push!(V2, x);
# Barsum cuemem x goals
# Barsum cuemem
X, (yval, ymin, ymax) = measure(R, "cuemem", stat=:mean, what=:amp, occNorm=false)
x =Gadfly.plot(X, scalecm, color=:cuemem, x=:cm, y=yval, ymin=ymin, ymax=ymax,
               Geom.bar, Geom.errorbar);
push!(V1, x);
# Barsum cuemem x goals
X, (yval, ymin, ymax) = measure(R, ["stopWell", "cuemem"], stat=:mean, what=:amp,
                    occNorm=false)
x = Gadfly.plot(X, x=:cm, color=:stopWell, y=:amp_mean, Geom.bar, Geom.errorbar);
push!(V2, x);
# Barsum cuemem
X, (yval, ymin, ymax) = measure(R, ["cuemem"], stat=:mean, what=:duration, occNorm=false)
x =Gadfly.plot(X, scalecm, color=:cuemem, x=:cm, y=yval, ymin=ymin, ymax=ymax,
               Guide.ylabel("duration_mean", orientation=:vertical), Geom.bar,
               Geom.errorbar)
push!(V1,x)
# TODO Barsum cuemem x goals :: Bootstraps have to be same elgnth
# X, (yval,ymin,ymax) = measure(R, ["startWell","stopWell","cuemem"], stat=:mean,
#                     what=:duration, occNorm=false)
# coord = Coord.cartesian(xmin=0.5, xmax=5.5, ymax=0.8)
# x = Gadfly.plot(X, scalecm, xgroup=:startWell, color=:cuemem, x=:stopWell,
#                y=yval, ymin=ymin, ymax=ymax, Geom.subplot_grid(coord,
#                                                                Geom.bar(position=:dodge),
#                                                                Geom.errorbar));
# push!(V3, x)
# Barsum cuemem
X, (yval,ymin,ymax) = measure(R, ["cuemem"], stat=:sum, what=:ampdur, occNorm=true)
x =Gadfly.plot(X, scalecm, color=:cuemem, x=:cm, y=yval, ymin=ymin, ymax=ymax,
               Guide.ylabel("Collective duration\n(occnorm)",
                            orientation=:vertical), Geom.bar, Geom.errorbar)
push!(V1,x)

V1,V2 = even_out(V1,V2)
set_default_plot_size(23cm,13cm)
p =vstack(hstack(V1...), hstack(V2...))

plotpath = plotsdir("ripples","summary_counts_and_amps_by_cuememStopWell_v2")
gsave(p, plotpath)


histogram(ripples.duration,label="Ripple duration")
savefig(plotsdir("ripples","duration.svg"))
