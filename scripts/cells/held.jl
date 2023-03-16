using GoalFetchAnalysis
using .Plot
using Plots, StatsBase, StatsPlots, DataFrames, DataFramesMeta, HypothesisTests

Plot.setparentfolder("cells")
Plot.setfolder("held")
Plot.setappend("animal=$(animal)_day=$(day)")
Plot.printstate()

# All epochs
# -----------
@df combine(groupby(combine(groupby(spikes, :unit), :epoch =>
                unique),:unit),nrow) begin
                histogram(:nrow, title="cells held for how many epochs")
end

# Run epochs
# -----------
@df combine(groupby(@subset(celleps, :task .== "cm"),:unit), nrow) begin
                histogram(:nrow, title="Held cells over runs")
end
heldrun = combine(groupby(@subset(celleps, :task .== "cm"),:unit), 
                nrow => :held_epochs)
Load.celltet.save_cell_taginfo(heldrun, "RY16", 36, "held_epochs")

# Prepare for a 10 minute bucketed approach
dT = diff(beh.time)
bins = Int(maximum(floor.(cumsum(dT)/60/10))) #get number of 10 minute buckets
beh.bins = DIutils.binning.digitize(beh.time, bins);
DIutils.filtreg.register(beh, spikes, on="time", 
    transfer=["bins"]);

# Get the number of spikes in each 10 minute bucket
unique_units= combine(groupby(spikes, :bins), :unit => x->length(unique(x)), renamecols=false)
@df unique_units plot(:bins, :unit, 
                title="Number of cells held for each 10 minute bucket")
Plot.save("held_cells per 10 minute bucket")

ep_unique_units= combine(groupby(spikes, :epoch), :unit => x->length(unique(x)), renamecols=false)
@df ep_unique_units plot(:epoch, :unit, 
                title="Number of cells held for each epoch")
vline!([ep_unique_units.epoch[argmax(ep_unique_units.unit)]])
Plot.save("held_cells per epoch")

# Get the firing rate of each cell in each 10 minute bucket
FR=  combine(groupby(spikes, [:bins, :unit]),
                :time => length => :spikecount,
                :time => first => :start,
                :time => last => :stop,
                renamecols=false)
transform!(FR, [:spikecount,:start,:stop] => ((x...) ->x[1]./(x[3].-x[2])) => :firingrate)
FR.firingrate =  replace(FR.firingrate, Inf=>missing)
dropmissing!(FR, :firingrate)
transform!(FR, :firingrate => (x->log10.(x)) => :logfiringrate)

# Unstack cells
FRun = unstack(FR, :bins, :unit, :firingrate)
logFRun = unstack(FR, :bins, :unit, :logfiringrate)

plot(
heatmap(Matrix(FRun[:, Not(:bins)]), title="Firing rate of held cells over time"),
heatmap(Matrix(logFRun[:, Not(:bins)]), title="LOG firing rate of held cells over time")
)
Plot.save("held_heatmap FR, 10 min bins")

