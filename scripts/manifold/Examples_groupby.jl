import Table
nExample = 20
groups   = [:startWell, :stopWell]



em = Table.convert_types.from_dimarrary(em, (beh, "time", ["startWell","stopWell","traj"]))
em = groupby(em, [:startWell,:stopWell, :traj])
plot(eachrow(hcat(em[3].data...))...)
