"""
Module for interacting with and applying
dynamic time warp to tabular data
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
    using Statistics, Interpolations, RollingFunctions
    using DataFrames
    using ProgressMeter
    using SearchSortedNearest

    export get_groupedexamples
    function get_groupedexamples(spikes, beh;
        dims = ["startWell","stopWell","traj"],
        values = ["x","y", "time"], thresh=12,
        printstats = true,
        groupingdim = :traj)::Matrix

        X = torate(spikes, beh)
        X = rate_todataframe(X, (beh,"time",[values...,dims...],))
        X = Table.group.equalize(X,  setdiff(dims,["traj"]), :traj; thresh)
        
        # Tensorize our dataframe
        tens = tensorize(X, dims, values)

        # And place the relevant grouping along the dimension
        grouped_examples = tenmat(tens, 
                findfirst(String.(dims) .== String.(groupingdim)))
        mis = vec((!).(any(ismissing.(grouped_examples), dims=1)))
        grouped_examples = grouped_examples[:,findall(mis)]
        grouped_examples = disallowmissing(grouped_examples)
        grouped_examples = [Matrix(g') for g in grouped_examples]
        
        if printstats
            nExamples, nGroups = size(grouped_examples);
            @info "stats" thresh nExamples nGroups
        end
        grouped_examples
    end

    export get_templates
    function get_templates(examples; method=:median)::Vector
        n,_ = size(examples[1])
        templates = []
        for i in 1:size(examples,2)
            data = examples[:,i][1:n,:]
            maxdatalen = maximum([size(d,2) for  d in data])
            data = [Matrix(d') for d in data]
            # Interpolate
            for i in 1:length(data)
                D = collect.(eachcol(data[i]))
                D = interpolate.(D, [BSpline(Linear())])
                D = hcat([d[LinRange(1,length(d), maxdatalen)]
                          for d in D]...)
                data[i]=D
            end
            template = median(cat(data...;dims=3),dims=3)[:,:]
            template = hcat([rollmedian(q,2) for q in eachcol(template)]...)
            push!(templates, template)
        end
        templates
    end

    export get_dtwtable
    function get_dtwtable(examples, templates)::DataFrame
        n,_ = size(examples[1])
        # APPLY TEMPLATES :: DO DTW
        warptablegroups = []
        @showprogress "dtw table for groups" for i in 1:size(examples,2)

            data = examples[:,i]
            inputdata = getindex.(examples[:,i], [1:n], [:])
            inputtemplates = templates[i][:,1:n]'
            dtwres = dtw.([inputtemplates], inputdata)
            seq1, seq2  = [a[2] for a in dtwres],
                          [a[3] for a in dtwres]
            
            template_time = templates[i][:,end]
            template_time .-= minimum(template_time)

            # Get time effects
            warptable = [DataFrame([template_time[s1] d[end, s2] s1 s2],["time_template_zerod","time","s1","s2"])
                         for (s1,s2,d) in zip(seq1, seq2, data)]
            transform!.(warptable, [:time => (x-> x .- minimum(x)) => :time_zerod])
            transform!.(warptable, [[:time_template_zerod,:time_zerod]=>((x,y)->(x .- y))=> :delta])
            warptable = vcat(warptable...; source=:warpexample)
            push!(warptablegroups, warptable)
        end
        W = vcat(warptablegroups...; source=:warpgroup)
        W[:, [:warpgroup, :warpexample, :time, :s1, :s2, :time_template_zerod,
              :time_zerod]]
    end

    export apply_warps
    function apply_warps(data::DataFrame, warptable::DataFrame;
                        on=:time, only_warp=true, nan_missing=true)::DataFrame
        if !issorted(data[!,on])
            sort!(data,on)
        end
        wgs_stats = if nan_missing
            combine(groupby(data,:warpgroup),
            :s1 => minimum, :s1 => maximum)
        else
            nothing
        end
        # Get closest match times (first for each index of the template)
        wgs = groupby(warptable, [:warpgroup, :warpexample])
        results=[]
        for wg in wgs
            push!(results, apply_warps_wg_first(data, wg))
        end
        vcat(results...)
    end
    function apply_warps_wg_first(data::DataFrame, wg::SubDataFrame;
                        nan_pad=nothing,
                        on=:time)::DataFrame

        wgf = combine(groupby(wg, :s1), first)
        wgf.time
        matches = searchsortednearest(data[!,on], wgf[!,on])
        data = data[matches,:]
        hcat(data, wgf[!,[:warpgroup, :warpexample, :time_template_zerod, :time_zerod]])
    end

    function warped_df_to_tensor(data::DataFrame, dims, vars; on=:time)
        # tensorize by dims and vars

    end

end

