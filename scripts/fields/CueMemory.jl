# ---------------------------------------------------------------------------
# CUE VERS MEMORY
# VEL > 1

memfields, memkdfields, grid = field.get_fields(beh, spikes; filters=filters)
cuefields, cuekdfields, grid = field.get_fields(beh, spikes; filters=filters)

# Generate all of the fields
# --------------------------
# Let's detect what's there
a=field.show_all_fields(cuekdfields, cells; textcolor=:black)
b=field.show_all_fields(memkdfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "kde_vel>1")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])
a=field.show_all_fields(cuefields, cells; textcolor=:black)
b=field.show_all_fields(memfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "hist_vel>1")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])

#                                              
# ,---.          |        ,---.o     |        |
# |  _.,---.,---.|        |__. .,---.|    ,---|
# |   ||   |,---||        |    ||---'|    |   |
# `---'`---'`---^`---'    `    ``---'`---'`---'
# ANY SPEED

filters = Dict("velVec"=> x->abs.(x).>=0,
               "cuemem"=> x->x.==1)
memfields, memkdfields = field.get_fields(beh, spikes;
                                               filters=filters,
                                               props=goalprops
                                              )
filters = Dict("velVec"=> x->abs.(x).>=0,
               "cuemem"=> x->x.==0)
cuefields, cuekdfields = field.get_fields(beh, spikes; 
                                               filters=filters,
                                               props=goalprops
                                              )

# Generate all of the fields
# --------------------------
# Let's detect what's there
a=field.show_all_fields(cuekdfields, cells; textcolor=:black)
b=field.show_all_fields(memkdfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "kde_vel>0")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])
a=field.show_all_fields(cuefields, cells; textcolor=:black)
b=field.show_all_fields(memfields, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "hist_vel>0")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])


