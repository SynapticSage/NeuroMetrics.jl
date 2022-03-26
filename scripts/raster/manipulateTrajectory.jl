"""
#+TITLE: manipulateTrajectory
#+PURPOSE: to quickly visualize the raster sorted by
#           behavioral properties
#
#+METHOD: Interact.jl - we will use this to quickly
#         replot the raster with behavioral overlays
#
# * Manipulandums
#   - BehavioralOverlay
#   - SortY = default()
#   - SortColor = default()
#   - ColorScheme = default()
#   - TimeGrouping = default()
#   - VBlock = default(traj)
#
# * Additional sortables need to add
#   - centroid for a list of behavioral properties
#     per cell
#
# * Plot funcs
#   - Raster for Gadfly.jl per time grouping
#
# * DIO additions to beh
#   - nose_in_well
#   - correct*nose_in_well
"""
module workspace
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Pushover, Revise
using Interact, Blink, Mux, ProgressMeter
using Statistics
__revise_mode__ = :eval
#includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("table.jl"))
includet(srcdir("raster.jl"))


animal = "RY16"
day=36
abbreviated = true
if abbreviated
    data_source=[ "spikes","behavior"]
    spikes, beh  = raw.load(animal, day; data_source=data_source)
else
    data_source=[ "spikes","behavior","lfp","cells","tetrode","ripples"]
    spikes, beh, lfp, cells, tetrode, ripples = raw.load(animal, day; 
                                                        data_source=data_source)
end
send(getPushoverClient(), "Finished loading julia data")


using Interact, Blink, Mux, ProgressMeter, Statistics, DataStructures
bprops = ["x", "y", "currentAngle", "currentPathLength", 
        "egoVec"];
splits= [nothing, "area", "ampEgoVec", "angEgoVec"]
split_type = Dict(nothing => :row,
                  "area"  =>:row,
                  "ampEgoVec" =>:column,
                  "angEgoVec" =>:column
                 )
timegroupings = ["traj", "subblock", "block"];
dio_props     = ["poke", "well"];
#sp = table.add_columns_from_othertable(spaugmented_propsikes, beh, bprops)
lookupcols = Dict((source=1,target=2) =>
                  table.expand_colnames(beh,union(bprops,timegroupings,dio_props)))
beh, spikes = raw.filterTables(beh, spikes; filters=Dict(),
                                  lookupcols=lookupcols)
spikes, sprops = table.handle_complex_columns(spikes, props_to_mod=bprops)
beh, xprops    = table.handle_complex_columns(beh, props_to_mod=bprops)
xprops = xprops[xprops.!="egoVec"]
spikes = raster.populate_sort_fields(spikes, beh, xprops);



## ------- ##
## BUTTONS ##
## ------- ##
val(x) = x.output.val;
behavior_vars   = val(dropdown(OrderedDict(zip(xprops, xprops)), 
                           value=xprops[1],
                        label="Behavior to plot"));
# Sort y axis
sortY = val(dropdown(OrderedDict(zip(xprops, xprops)), 
                     value=xprops[1],
                        label="Sort Y by =>"));
# Sort color axis
sortColor = val(dropdown(OrderedDict(zip(xprops, xprops)), 
                         value=xprops[1],
                        label="Color spikes by =>"));
# Split data by
subplot  = val(dropdown(OrderedDict(zip(splits, splits)), 
                    value="area",
                    label="Subplot by =>"));
# Name of time grouping
timegroup = val(togglebuttons(timegroupings,
                             value="traj",
                             label="Color time group=>"));
# Name of time grouping
colorTimeGroup = val(togglebuttons(timegroupings,
                             value="traj",
                             label="Color time group=>"));
# Groups of time
uTraj = sort(unique(spikes[:, timegroup]));
timegroup_count = val(dropdown(OrderedDict(zip(string.(uTraj), string.(uTraj))),
                        label="Trajectory",
                        value="1"));
toggle_subplot_rank   = val(togglebuttons(["none", "sortby", "colorby", "both"],
                              value="none",
                              label="Use area-specific rank order"));
# Metric/Rank toggle
toggle_metricRank = val(togglebuttons(["none", "sortby", "colorby", "both"],
                              value="none",
                              label="Uses metric instead of rank order"));

# Other testing variables
rankOrder, rankColor, rankY = toggle_metricRank, nothing, nothing

## ------- ##
## TEST    ##
## ------- ##
w = Window();
ui = @manipulate for timegroup in timegroups, timegroup_count in timegroup_counts, sortY in sortYs, sortColor in sortColors, subplot in subplots, srToggle in toggle_subplot_rank in mrToggle in toggle_metricRank, behavior_vars in behavior_vars

raster.gadfly_subplot(spikes, beh, 
       timegroup, parse.(Int,timegroup_count); # TIME CHUNKS
       sortY=sortY, # Sort Y axis how?
       sortColor=sortColor,  # Color points by
       subplot=subplot, # Create subplots by
       rankOrder=rankOrder,
       rankColor=rankColor,
       rankY=rankY, # EVERYTHING AFTER THIS INTO (OTHER...) construct
       group_tabledim=:row,
       stackdir="v",
       colorTimeGroup=colorTimeGroup, # Geom.HBand groups of time
       behavior_vars=behavior_vars, # Geom.line behavior over raster
       );

end
body!(w,ui)
