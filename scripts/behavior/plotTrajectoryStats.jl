includet(srcdir("utils.jl"))
import .utils
import Plots

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

# 
