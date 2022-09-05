@time spikes, beh, ripples, cells = Load.load("RY16", 36);
Load.filtreg.register(beh,spikes,on="time", transfer=["traj"])

Dict( g.traj[1] => (p=plot(g.t, g.x); plot!(g.t,g.y); p)
     for g in groupby(beh, :traj))

