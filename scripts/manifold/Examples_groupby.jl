import Table
em = embedding[key]
nExample = 20
groups   = [:startWell, :stopWell]
em = Table.from_dimarray(em, (beh, "time", String.(groups)))
