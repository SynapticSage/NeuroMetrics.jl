# ---------------------------------------------------------------------------
# CUE VERS MEMORY
# VEL > 1
mem = field.get_fields(beh, spikes; filters=filters)
cue = field.get_fields(beh, spikes; filters=filters)

# Generate all of the fields
# --------------------------
# Let's detect what's there
a=field.show_all_fields(cue.kR, cells; textcolor=:black)
b=field.show_all_fields(mem.kR, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "kde_vel>1")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])
a=field.show_all_fields(cue.kR, cells; textcolor=:black)
b=field.show_all_fields(mem.kR, cells; textcolor=:black)
overall = Plots.plot(a,b)
plots = plotsdir.("fields", "goalfield" .* ["cue", "mem", "cuemem"] .* "hist_vel>1")
sf(a, plots[1])
sf(b, plots[2])
sf(overall, plots[3])


