
# All epochs
# -----------
@df combine(groupby(combine(groupby(spikes, :unit), :epoch => unique),:unit),nrow) histogram(:nrow, title="cells held for how many epochs")

# Run epochs
# -----------
@df combine(groupby(@subset(celleps, :task .== "cm"),:unit), nrow) histogram(:nrow, title="Held cells over runs")
heldrun = combine(groupby(@subset(celleps, :task .== "cm"),:unit), nrow => :held_epochs)
Load.celltet.save_cell_taginfo(heldrun, "RY16", 36, "held_epochs")

