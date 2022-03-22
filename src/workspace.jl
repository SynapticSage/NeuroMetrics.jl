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

function getPushoverClient()

    token = open(expanduser("~/.pushover.token"),"r") do f
        token = read(f, String)
    end
    user = open(expanduser("~/.pushover.user"),"r") do f
        user = read(f, String)
    end

    return PushoverClient(user, token)

end

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



end

