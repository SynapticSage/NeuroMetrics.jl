@time using GoalFetchAnalysis
using Test, Revise
using DrWatson

watched = keys(Revise.watched_files)
#dirs = ["", "Utils", "Table", "Load", "Field", "Timeshift", "Munge", "Decode"]
dirs = ["", "Utils", "Table", "Load", "Field", "Timeshift", "Munge"]
dirs = srcdir.(dirs)

# Most basic test! Do my libraries wind up in Reivse's watched files?
# If not, you may have either a compile fail or a flaw in the DAG-like
# dependency structure
for dir in dirs
    if endswith(dir, '/')
        dir = dir[1:end-1]
    end
    dir_in_watched = dir ∈ watched
    !(dir_in_watched) ? @warn("dir=$dir not in watched") : nothing
    @test dir_in_watched
end


#     --.--,---.,---.--.--    ,---.,---.. . .   #
#       |  |--- `---.  |      |---'|---|| | |   #
#       |  |        |  |      |  \ |   || | |   #
#       `  `---'`---'  `      `   ``   '`-'-'   #
#println("Loading data")
#spikes   = raw.load_spikes("RY16", 36)
#behavior = raw.load_behavior("RY16", 36)
#
#println("Testing registration")
#register_prop = ["velVec", "x", "y"];
## Test registrattion
#behavior, spikes = raw.register(behavior, spikes; 
#                                transfer=register_prop, 
#                                on="time");
#@test all(in.(register_prop, [names(spikes)]))
#
## Test filtration
#filters = Dict("velVec"=> vel->abs.(vel) .> 5,
#     "x"     => x->(x.>2) .& (x .< 400))
#sp = raw.filter(spikes, filters=filters)[1]
#@test all(abs.(sp.velVec) .> 5)
#@test all((sp.x .> 2) .& (sp.x .< 400))
#
## Test mergence of the two
#println("raw.filterTables")
#spikes, behavior = raw.filterTables(spikes, behavior; filters=filters, lookupon="time", lookupcols=register_prop)
#@test all(abs.(spikes.velVec) .> 5)
#@test all((spikes.x .> 2) .& (spikes.x .< 400))
#println("raw.filterAndRegister")
#spikes, behavior = raw.filterAndRegister(spikes, behavior; filters=filters, on="time", transfer=register_prop)
#@test all(abs.(spikes.velVec) .> 5)
#@test all((spikes.x .> 2) .& (spikes.x .< 400))
#
## TODO Test transfer can be given as a Dict{NamedTuple, Vector{String}} instruction
#
#
#
## --.--          |        ,---.o     |        |
##   |  ,---.,---.|---     |__. .,---.|    ,---|
##   |  |---'`---.|        |    ||---'|    |   |
##   `  `---'`---'`---'    `    ``---'`---'`---'
#props = ["x", "y"]
#kws=(;resolution=80, splitby=[:unit,:area], filters=merge(filt.speed_lib, filt.cellcount))
#newkws = (; kws..., gaussian=2.3*0.5, props=props, filters=merge(kws.filters))
#place = field.get_fields(behavior, spikes; newkws...);
