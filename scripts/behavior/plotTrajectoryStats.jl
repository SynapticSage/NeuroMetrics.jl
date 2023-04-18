includet(srcdir("utils.jl"))
import .utils
import Plots, DataFrames, DataFramesMeta
using DIutils
using DI
using ImageFiltering
using GoalFetchAnalysis
import GoalFetchAnalysis.Plot

opt["animal"] = "animal" in keys(opt) ? opt["animal"] : "RY16"
opt["day"]    = "day"    in keys(opt) ? opt["day"]    : 36

Plot.setparentfolder("behavior")
Plot.setfolder("trajectory")
Plot.setappend("animal=$(opt["animal"])_day=$(opt["day"])")

spikes, beh = DI.load(opt["animal"], opt["day"],
                     data_source=["spikes", "behavior"])
animal, day = opt["animal"], opt["day"]

G = combine(groupby(beh, :traj),
          :time=> (x->diff([extrema(x)...])) => :range).range;
Plots.histogram(G, bins=20000, xlim=(0,1), label="traj time dur", xlabel="min");
vline!([median(G)], c=:black, linestyle=:dash, linewidth=3, label="typical");
annotate!(median(G)*2, 20, "typical\n$(round(median(G),digits=3))" )

savefig(plotsdir("behavior","trajectory","typical_trajectory_length.png"))
savefig(plotsdir("behavior","trajectory","typical_trajectory_length.svg"))
savefig(plotsdir("behavior","trajectory","typical_trajectory_length.pdf"))

task    = mean(   (!).((isnan.(beh.traj)))    )
drinkingAndNonTask = 1-task
p1 = Plots.pie(["Task", "Drinking+NonTask"], [task, drinkingAndNonTask])

task    = mean(  (!).(isnan.(beh.traj)) .&& (abs.(beh.velVec) .> 2) )
drinkingAndNonTask = 1-mean(isnan.(beh.traj))
p2 = Plots.pie(["Running task (>2cm/s)", ""], [task, drinkingAndNonTask])

task    = mean(  (!).(isnan.(beh.traj)) .&& (abs.(beh.velVec) .> 5) )
drinkingAndNonTask = 1-mean(isnan.(beh.traj))
p3 = Plots.pie(["Running task (>5cm/s)",""], [task, drinkingAndNonTask])

task    = mean(  (!).(isnan.(beh.currentHeadEgoAngle)) )
drinkingAndNonTask = 1-mean(isnan.(beh.traj))
p4 = Plots.pie(["Goal-path-angle times",""], [task, drinkingAndNonTask])

Plots.plot(p1, p2, p3, p4)
savefig(plotsdir("behavior","trajectory","typical_trajectory_piecharts.png"))
savefig(plotsdir("behavior","trajectory","typical_trajectory_piecharts.svg"))
savefig(plotsdir("behavior","trajectory","typical_trajectory_piecharts.pdf"))

# RELATIVE TRAJECTORY TIMES
Plots.histogram(beh.trajreltime[beh.stopWell.!=-1], label="Relative trajectory
                time", xlabel="Relative traj\n0:start 1:stop",
                ylabel="samples")
utils.savef("behavior","trajectory","histogram_relative_trajectory_times")

beh = groupby(beh, :stopWell)
edges = 0:0.1:1
centers = mean([edges[1:end-1] edges[2:end]], dims=2)
p = []
for b in beh
    i = b.stopWell[1]
    heights = fit(Histogram, b.trajreltime, edges).weights
    (m,M) = extrema(heights)
    push!(p,
          Plots.bar(centers, heights, label="$i", ylim=(m,M),
                            xlabel="Relative traj\n0:start 1:stop",
                            ylabel="samples"))
end
beh = combine(beh, identity)
Plots.plot(p...)
utils.savef("behavior","trajectory","histogram_by_stopWell_relative_trajectory_times")

# DOes this evolve over trajectory
Plots.histogram2d(beh.trajreltime, beh.traj, label="Relative trajectory time",
                  xlabel="Relative traj\n0:start 1:stop", ylabel="traj")
utils.savef("behavior","trajectory","histogram2d_relative_trajectory_times_versus_traj")


# Typical time length when behavior > 2cms
histogram(trajperiod.δ .* trajperiod.frac, xlim=(0,20))
xlims!(0,25)
m=median(trajperiod.δ .* trajperiod.frac)
M=mean(trajperiod.δ .* trajperiod.frac)
vline!([median(trajperiod.δ .* trajperiod.frac)], xlim=(0,20), c=:black, linestyle=:dash, title="Traj lengths in time\nmedian=$m\nmean=$M")
savefig(plotsdir("behavior","trajectory-timing","times,gt2cms.png"))
savefig(plotsdir("behavior","trajectory-timing","times,gt2cms.pdf"))
# Have to throw out the outliers if want a realistic mean. Median works.

# ------------------------------------
# PLOT: Typical trajectory speed curves
# ------------------------------------
# Speed changes
diff(beh.speed) |> histogram
using Loess, RollingFunctions
# l = loess(beh.time, beh.speed, span=0.1)
# beh[!,:speedsmooth] = predict(l, beh.time)
k=12
beh.speedsmooth = [zeros(Int(round(k/2))-1); 
    rollmean(beh.speed, k); zeros(Int(k/2))]
begin 
    tmp=@subset(beh, :traj .== 2) 
    @df tmp plot(:time, :speed)
    @df tmp plot!(:time, :speedsmooth)
    plot!(xlabel="Time (s)", ylabel="Speed (cm/s)")
    ylims!(0,10)
end
beh.movingsmooth = beh.speedsmooth .> 2.5
beh.moving  = beh.movingsmooth

# ------------------------------------
# PLOT: Typical trajectory speed curves
# ------------------------------------
b   = Base.filter(x->x.traj != NaN && 
                  x.traj !== missing, beh)
color = DIutils.plotutils.getplotcolor(b.blocktraj, :vik)
hatch = map(b.correct) do x
    if x == 1
        return nothing
    elseif x == 0
        return :\
    else
        return :x
    end
end
b.color = color
b.hatch = hatch
b = groupby(b, :traj)
P = [plot(b[i].time, abs.(b[i].speedsmooth), 
    c=b[i].color, fillstyle=b[i].hatch, fill=0,
    xticks=(extrema(b[i].time)|>collect, 
            round.(extrema(b[i].time), digits=2)|>collect),
) for i in eachindex(b)]

layout = @layout [a{0.005h}; grid(7,7)]
blank = Plot.blank(title="Typical trajectory speed curves"*
    "animal=$animal, day=$day",
    titlefontsize=8, legend=:none);

C = []
for (i,c) in enumerate(2:(7*7):length(P))
    m = min(c+(7^2-1), length(P))
    kws = i==1 ? (;layout) : (;)
    p=plot(blank, P[c:m]..., tickfontsize=3, label="", 
        margin=0.001Plots.mm; kws...)
    Plot.save("trajectory_collection=$i")
    push!(C, p)
end

# ---------------------------------------------------------------
# PLOT: Heatmap counting trajectories of each type (startWell, stopWell)
# ---------------------------------------------------------------
Plot.setfolder("start->stopwell")
# How many trajectories of each type?
epoch = groupby(beh, :epoch)
P = [] 
e = epoch |> first
for e in epoch
    b = @subset(e, 
                :traj .!= NaN,
                :startWell .!= -1, :stopWell .!= -1; view=true)
    hw = b.homewell[1]
    b = groupby(b, [:startWell, :stopWell])
    b = combine(b, :traj => (t->length(unique(t))) => :ntraj,
                  :time => (t->length(t)) => :ntime)
    b = sort(b, [:startWell, :stopWell])
    btraj = unstack(b, :startWell, :stopWell, :ntraj)
    btraj = btraj[!, sort(names(btraj))]
    btime = unstack(b, :startWell, :stopWell, :ntime)
    btime = btime[!, sort(names(btime))]
    p=Plots.heatmap(Matrix(btraj[:,Not(:startWell)]), 
        title="epoch=$(e.epoch[1])", 
        xlabel="Stop well", ylabel="Start well", aspectratio=1,
        xticks=1:5, yticks=1:5, xtickfont=Plots.font(8), ytickfont=Plots.font(8),
    xlim=(0.5,5.5), ylim=(0.5,5.5), c=:viridis)
    Plots.annotate!(1:5, fill(hw, 5), Plots.text("*", :black, :center))
    Plots.annotate!(fill(hw, 5), 1:5, Plots.text("*", :black, :center))
    Plot.save("heatmap_ntraj,epoch=$(e.epoch[1]),homewell=$hw")
    push!(P, p)
end
blank = Plot.blank(title="Heatmap counting trajectories of each type\n(startWell, stopWell)\n"*
    "animal=$animal, day=$day",
    titlefontsize=8, legend=:none);
layout = Plots.@layout [a{0.005h}; Plots.grid(1,length(epoch))]
Plots.plot(blank, P...; layout, textfontsize=3, tickfontsize=3, label="", 
    margin=0.001Plots.mm) 
Plot.save("heatmap_ntraj,all_epochs")

# ---------------------------------------------------------------
# PLOT: Heatmap counting trajectories of each type (startWell, stopWell)
# ---------------------------------------------------------------
# Home-arena splits per epoch
# THse suggest that the trials are innapropriately labeled!!! 
# ---------------------------------------------------------------
epoch = groupby(dropmissing(beh,:ha), [:epoch, :ha])
P = [] 
e = epoch |> first
for e in epoch
    b = @subset(e, :traj .!= NaN, :startWell .!= -1, :stopWell .!= -1;
                    view=true)
    hw = b.homewell[1]
    b = groupby(b, [:startWell, :stopWell])
    b = combine(b, :traj => (t->length(unique(t))) => :ntraj,
                  :time => (t->length(t)) => :ntime)
    btraj = unstack(b, :startWell, :stopWell, :ntraj)
    btraj = btraj[!, sort(names(btraj))]
    btime = unstack(b, :startWell, :stopWell, :ntime)
    btime = btime[!, sort(names(btime))]
    xticks = names(btraj[:,Not(:startWell)])
    xticks = (1:length(xticks), xticks)
    yticks = btraj.startWell
    yticks = (1:length(yticks), yticks)
    p=nothing
    try
        p=Plots.heatmap(Matrix(btraj[:,Not(:startWell)]), 
            title="epoch=$(e.epoch[1])\nha=$(e.ha[1])", 
            xlabel="Stop well", ylabel="Start well", aspectratio=1,
            xticks=xticks, yticks=yticks, xtickfont=font(8), ytickfont=font(8),
            xlim=(0.5,5.5), ylim=(0.5,5.5), c=:viridis)
        xhwloc = findall(parse.(Int,xticks[2]) .== hw)
        if !isempty(xhwloc)
            annotate!(fill(xhwloc[1], length(xticks[1])), yticks[1], Plots.text("*", :black, :center))
        end
        yhwloc = findall(yticks[2] .== hw)
        if !isempty(yhwloc)
            annotate!(xticks[1], fill(yhwloc[1], length(xticks[1])), Plots.text("*", :black, :center))
        end
        Plot.save("ha=$(e.ha[1])_heatmap_ntraj,epoch=$(e.epoch[1]),homewell=$hw")
    catch
        @infiltrate
    end
    push!(P, p)
println("length(P): ", length(P))
end
blank = Plot.blank(title="Heatmap counting trajectories of each type\n(startWell, stopWell)\n"*
    "animal=$animal, day=$day",
    titlefontsize=8, legend=:none);
layout = @layout [a{0.005h}; grid(2,Int(length(epoch)/2))]
plot(P...; textfontsize=3, tickfontsize=3, label="", 
    margin=0.001Plots.mm)

# ---------------------------------------------------------------
# PLOT: Heatmap of number of trajectories per epoch
# ---------------------------------------------------------------
# Let's extract immobility periods ...
# For every trajectory, there will be a single mobility period and
# a single immobility period.
# ---------------------------------------------------------------
DI.smooth_movement!(beh)
beh[!,:movingsmooth] = imfilter(beh.moving_speedsmooth,
    Kernel.gaussian((movebool_gaussian*10,)))
B = groupby(
    @subset(beh, :traj .!= NaN, :startWell .!= -1, :stopWell .!= -1; view=true),
    :traj)
P = []
for b in B
    p = 
    plot(
        (plot(b.speedsmooth, label="smoothed", fill=0, fillalpha=0.5);
         hline!([4], label="threshold", color=:black, linestyle=:dash);
            plot!(b.speed,label="speed")),
        (plot(b.movingsmooth,label="smoothed", fill=0, fillalpha=0.5);
         hline!([0.5], label="threshold", color=:black, linestyle=:dash);
         plot!(b.moving,label="moving"));
        layout=grid(2,1), legendposition=:outerbottomright,
        title="example of how handling mobility"
    )
    push!(P, p)
end

beh[!,:immobility] = beh.movingsmooth .< 0.5
b = @subset(beh, :traj .!== missing,
    :traj .!= NaN, :startWell .!= -1, :stopWell .!= -1; view=true
)
b = groupby(b, :traj)

# ---------------------------------------------------------------
# PLOT: speed, immobility, delta_xy
# ---------------------------------------------------------------
P=[]
t = b |> first
delta_xy(x1,x2) = sqrt(sum((x2-x1).^2))
for t in b
    p=plot(t.time, t.speed, label="speed")
    plot!(t.time, t.immobility .* ylims()[2], fill=0, label="immobility",
        fillalpha=0.2, 
    )
    plot!(t.time, delta_xy.(eachrow([t[1,:x] t[1,:y]]), eachrow([t.x t.y])),
        fill=0, c=:green, fillalpha=0.2,
        label="delta_xy")
    push!(P, p)
end
plot(P[1:24]..., layout=grid(3,8), size=(2000,1000), legendbackgroundalpha=0.5)


