module workspace
using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using Revise
using Interact, Blink, Mux, ProgressMeter
using Statistics
__revise_mode__ = :eval
#includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("table.jl"))
includet(srcdir("raster.jl"))
includet(srcdir("utils.jl"))



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

send(utils.getPushoverClient(), "Finished loading julia data")



end

