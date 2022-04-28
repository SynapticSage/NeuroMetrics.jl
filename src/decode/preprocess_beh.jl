export annotate_pastFutureGoals

function annotate_pastFutureGoals(beh)

    beh = groupby(beh, [:epoch, :traj])
    for (g,group) in enumerate(beh)
        if g!=length(beh)
            group.futureStopWell .= beh[g+1].stopWell[1]
        end
        if g!=1
            group.pastStopWell .= beh[g-1].stopWell[1]
        end
    end
    beh = sort(combine(beh, identity), :time)
    replace!(beh.pastStopWell, missing=>-1)
    replace!(beh.futureStopWell, missing=>-1)
    if doPrevPast
        beh.hw = argmax(StatsBase.fit(StatsBase.Histogram,
                                      filter(b->b!=-1,beh.stopWell), 1:6).weights)
        beh = groupby(beh, [:epoch, :block])
        for (g,group) in enumerate(beh)
            # TODO
        end
        replace!(beh.prevPastStopWell, missing=>-1)
        beh = sort(combine(beh, identity),:time)
    end
end

function annotate_poke(beh)
    pn = sort([name for name in names(beh) if occursin("poke_", name)])
    poke_matrix = replace(Matrix(beh[!, pn]), NaN=>0)
    poke_matrix = BitMatrix(poke_matrix)
    pn = replace([findfirst(row) for row in eachrow(poke_matrix)],nothing=>0)
    beh.poke = pn
    beh
end
