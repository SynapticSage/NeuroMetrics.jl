if isdefined(Main, :DrWatson)
    include(scriptsdir("reactivate","imports.jl"))
else
    include("imports.jl") # include(scriptsdir("reactivate","imports.jl"))
end

# Get and prep the data!
# ------------------------------------------------------------------
spikes = DI.load_spikes(opt["animal"],   opt["day"])
beh    = DI.load_behavior(opt["animal"], opt["day"])
cells  = DI.load_cells(opt["animal"],    opt["day"])
# Get firing rate DimArray
R = spiking.torate(spikes, beh, gaussian=0.20)
# Create dataframe of R
Rdf = DataFrame(R, name=:rate)
DIutils.pushover("Finished creating Rdf")
# Create a new field enumerating (startWell, stopWell) combinations
# with integers
S  = eachrow(Matrix(beh[:, [:startWell, :stopWell]])) |> collect
uS = unique(S, dims=1)
sf = [findfirst(isequal(s), uS) for s in S]
beh[:, :startstopWell] = sf
# Create a field specifying whether the animal is still or moving
beh[:, :moving] = beh[:, :speed] .> 2 # cm/s
groupings = [:startstopWell, :startWell, :stopWell, :moving, :ha, :correct, 
    :epoch]
println("Available trajectories:", unique(beh[:, :traj]))

# Register that to Rdf
sort!(Rdf, [:time, :unit])
register(beh, Rdf,   on="time", transfer=["traj",string.(groupings)...])
Rdf = @subset(Rdf, :stopWell .!= -1, :startWell .!= 1)
sort!(Rdf, [:unit, :time])
register(cells, Rdf, on="unit", transfer=["area"])
sort!(Rdf, [:time, :unit])
ca1, pfc = groupby(Rdf, :area)[1:2];
@assert(ca1.area[1] == "CA1" && pfc.area[1] == "PFC", 
"Areas are not CA1 and PFC")
ca1 = ca1[:, Not(:area)];
pfc = pfc[:, Not(:area)];
Rdf = Rdf[:, Not(:area)];
dropmissing!(ca1, :startstopWell)
dropmissing!(pfc, :startstopWell)
DIutils.pushover("Finished creating area splits")
ca1.ha = replace(ca1.ha, missing => 'm')
ca1    = groupby(ca1, groupings);
pfc    = groupby(pfc, groupings);
# INFO: up to this point, all startWell elements exist!
println("Available trajectories, Rdf:", unique(Rdf[:, :traj]))
println("Available trajectories, ca1:", unique(combine(ca1,identity)[:, :traj]))

# Speed changes
diff(beh.speed) |> histogram

# Doubounce movement
DI.debounce_movement!(beh; threshold=0.2)
beh.moving = beh.movingdeb

# IF startWell=1 is actually, being rejected, I need to actually plot the conditions
# that lead to its rejection
# @subset(Rdf, :startWell .== 1) |> Voyager()
# @subset(beh, :startWell .== 1) |> Voyager()

GC.gc()

# ------------------------------------------------------------------
# Which data will we accept?, right now I'm only accepting combinations
# that have more than num_cells samples
# ------------------------------------------------------------------
# BUG:
# why so many 'h' 'a' missing? especially when start and stopwell are both
# defined? 
# It turns out ALL of the missing values hail from the incorrect trials
# ------------------------------------------------------------------
accepted = begin
    # This line works as follows: 
    # 1. Take the first row of each cell group
    # 2. Filter out the groups that have less than num_cells samples
    # 3. Filter out the groups that have less than 3 trajectories
    # 4. Map the groups to a new dataframe with the groupings as columns
    #    and the size of the group and number of unique trajectories in the
    #    group as new columns
    accepted = 
        Base.map(
        Base.filter(map(x->groupby(x,:unit)|>first ,collect(ca1))) do x #1
        size(x,1) > size(@subset(cells,:area.=="CA1"),1) && length(unique(x.traj)) ≥ 3 #2,3
    end) do x #4
        sz = size(x,1)
        ntraj = length(unique(x.traj))
        x=DataFrame(x[1,groupings])
        cols=propertynames(x)
        x[!,:size] .= sz
        x[!,:ntraj] .= ntraj
        x[!, [:size, :ntraj, cols...]]
    end
    accepted = sort(vcat(accepted...), :size, rev=true)
    accepted = dropmissing(accepted)
    # accepted |> Voyager()
        p=@df accepted scatter(:correct, :ha)
        p.subplots[1].series_list[1][:x] .+= 
                    0.02randn(size(p.subplots[1].series_list[1][:x]))
        p.subplots[1].series_list[1][:y] .+= 
                    0.05randn(size(p.subplots[1].series_list[1][:y]))
        plot(p)
    unicodeplots()
        histogram((accepted.size * 0.033333)./60)
        histogram(accepted.ntraj, xlims=(0,20))
        histogram(combine(groupby(accepted,[:startstopWell,:epoch]),
                :ntraj=>sum).ntraj_sum)
    gr()
    accepted
end
# INFO: `accepted` :: there are no accepted startWell=1 combinations for RY16,36


# ------------------------------------------
# Recombine ca1 and pfc, filter, and re-split
# ------------------------------------------
@time begin
    ca1 = combine(ca1, identity);
    pfc = combine(pfc, identity);
    inds = collect(eachrow(Matrix(ca1[!,groupings]))) .∈ 
          (collect(eachrow(Matrix(accepted[!,groupings]))),);
    inds = @time disallowmissing(replace(inds, missing=>false));
    print("Fraction of points kept: ", sum(inds)/length(inds))
    ca1  = ca1[inds,:];
    inds = collect(eachrow(Matrix(pfc[!,groupings]))) .∈ 
          (collect(eachrow(Matrix(accepted[!,groupings]))),);
    inds = disallowmissing(replace(inds, missing=>false))
    print("Fraction of points kept: ", sum(inds)/length(inds))
    pfc = pfc[inds,:];
    ca1 = groupby(ca1, groupings);
    pfc = groupby(pfc, groupings);
    if size(ca1,1) != size(pfc,1)
        println("size(ca1,1)=$(size(ca1,1)) != size(pfc,1)= $(size(pfc,1))")
        println("...Aligning!")
        K1, K2 = keys(ca1), keys(pfc)
        K = intersect(K1, K2)
        println("...Keeping $(length(K)) keys")
        println("---------------")
        println(K)
        println("---------------")
        pfc′, ca1′ = Dict(), Dict()
        for k in K
            k = NamedTuple(k)
            pfc′[k] = pfc[k]
            ca1′[k] = ca1[k]
        end
        ca1 = vcat(values(ca1′)...)
        pfc = vcat(values(pfc′)...)
    end
    ca1, pfc = groupby(ca1, groupings), groupby(pfc, groupings);
    @assert(size(ca1,1) == size(pfc,1),
    "size(ca1,1)=$(size(ca1,1)) != size(pfc,1)=$(size(pfc,1))")
    # @assert(map(
    # (x,y) -> size(x,1) == size(y,1), collect(ca1), collect(pfc)),
    # )
end
println("Available trajectories, ca1:", 
    unique(combine(ca1,identity)[:, :traj]))
println("Number of trajectories, ca1:", 
    length(unique(combine(ca1,identity)[:, :traj])))

include("react.jl")
