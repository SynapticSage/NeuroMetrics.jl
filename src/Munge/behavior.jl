module behavior

    export annotate_pastFutureGoals
    export annotate_relative_xtime!
    using StatsBase
    using DataFrames
    using ProgressMeter
    using Missings
    using Infiltrator
    using DataStructures: OrderedDict

    import DIutils, DI

    function annotate_pastFutureGoals(beh::DataFrame; doPrevPast::Bool=false)

        beh = groupby(beh, [:epoch, :traj])
        for (g,group) in enumerate(beh)
            if g!=length(beh)
                group.futureStopWell .= beh[g+1].stopWell[1]
                group.futureCuemem   .= beh[g+1].cuemem[1]
            end
            if g!=1
                group.pastStopWell .= beh[g-1].stopWell[1]
                group.pastCuemem .= beh[g-1].cuemem[1]
            end
        end
        beh = sort(combine(beh, identity), :time)
        replace!(beh.pastStopWell, missing=>-1)
        replace!(beh.futureStopWell, missing=>-1)
        if doPrevPast
            beh.hw .= argmax(StatsBase.fit(StatsBase.Histogram,
                                          filter(b->b!=-1,beh.stopWell), 1:6).weights)
            beh = groupby(beh, [:epoch, :block])
            for (g,group) in enumerate(beh)
                # TODO
            end
            replace!(beh.prevPastStopWell, missing=>-1)
            beh = sort(combine(beh, identity),:time)
        end
        beh
    end

    function annotate_poke!(beh::DataFrame)::DataFrame
        pn = sort([name for name in names(beh) if occursin("poke_", name)])
        poke_matrix = replace(Matrix(beh[!, pn]), NaN=>0)
        poke_matrix = BitMatrix(poke_matrix)
        pn = replace([findfirst(row) for row in eachrow(poke_matrix)],nothing=>0)
        beh[!,:poke] = pn
        beh
    end

    function annotate_relative_xtime!(beh::DataFrame, x=:traj, on=:time)::DataFrame
        find_relative(t,m,M) = (t.-m)./(M-m)
        beh = groupby(beh, x)
        @showprogress 0.1 "adding rel $x on $on" for group in beh
            m, M = extrema(group[!,on])
            group[!, String(x)*"rel"*String(on)] = find_relative(group[!,on], m, M)
        end
        combine(beh, identity)
    end


    export get_wells_df
    function get_wells_df(animal,day,epoch;task=nothing,beh=nothing)
        task = task === nothing ? DI.load_task(animal, day) : task
        beh = beh === nothing ? DI.load_behavior(animal, day) : beh
        beh = dropmissing(beh, :stopWell)
        wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
        wellrange = maximum(skipmissing(beh.stopWell))
        wells.well = collect(1:wellrange)
        wells.name .= ""
        homeWell, arenaWells = begin
            hw = argmax(StatsBase.fit(StatsBase.Histogram, filter(b->b!=-1,beh.stopWell), 1:wellrange).weights)
            @info "Homewell = $hw"
            wells[hw,:], wells[setdiff(1:wellrange,hw),:]
        end
        homeWell.name = "home"
        arenaWells.name .= "arena"
        wells = vcat(DataFrame(homeWell),arenaWells)
    end


end
