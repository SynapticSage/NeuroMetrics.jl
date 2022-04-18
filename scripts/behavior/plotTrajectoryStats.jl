
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
