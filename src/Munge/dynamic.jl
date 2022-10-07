"""
Module for interacting with dynamic time warp
of data
"""
module dynamic

    using DynamicAxisWarping
    using GoalFetchAnalysis
    using Munge.spiking
    using Table
    using Missings
    using Plot.task
    using Infiltrator
    using RecipesBase
    using TensorToolbox
    using Munge.tensor

    function get_groupedexamples()::Matrix
        unique(beh.subblock)
        unique(beh.traj)
        dims = ["startWell","stopWell","traj"]
        values = ["x","y", "time"]
        groupingdim = :traj
        n = length(values)-1
        X = torate(spikes, beh)
        X = rate_todataframe(X, (beh,"time",[values...,dims...],))
        X = Table.group.equalize(X,  setdiff(dims,["traj"]), :traj, thresh=12)
        
        # Tensorize our dataframe
        out = tensorize(X, dims, values)

        # And place the relevant grouping along the dimension
        out = tenmat(out, findfirst(String.(dims) .== String.(groupingdim)))
        mis = vec((!).(any(ismissing.(out), dims=1)))
        out = out[:,findall(mis)]
        out = disallowmissing(out)
        out = [Matrix(o') for o in out]
        
        printstats = true
        if printstats
            nExamples, nGroups = size(out);
            @info "stats" thresh nExamples nGroups
        end
        out
    end

    function get_templates()::Vector
    end

    function get_dtwtable()::DataFrame


end

