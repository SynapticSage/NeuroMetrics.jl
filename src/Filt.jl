module Filt

    using DataFrames
    using DataStructures
    export get_filters
    using Base: merge
    using Infiltrator

    function SPEED(x)
        abs.(x) .> 4
    end
    function SPEED_LIB(x)
        abs.(x) .> 2
    end
    function STILL(x)
        abs.(x) .> 2
    end
    speed     = OrderedDict("velVec"=>SPEED)
    speed_lib = OrderedDict("velVec"=>SPEED_LIB)
    still     = OrderedDict("velVec"=>STILL)

    function CORRECT(x)
        x.==1
    end
    function INCORRECT(x)
        x.==0
    end
    function NONTASK(x)
        (x.!=0) .&& (x.!=1)
    end
    function TASK(x)
        (x.==0) .|| (x.==1)
    end
    correct   = OrderedDict("correct" => CORRECT)
    incorrect = OrderedDict("correct" => INCORRECT)
    nontask   = OrderedDict("correct" => NONTASK)
    task      = OrderedDict("correct" => TASK)
    # Alias
    error     = incorrect

    function MEM(x)
        x.==1
    end
    function CUE(x)
        x.==0
    end
    cue       = OrderedDict("cuemem" => CUE)
    mem       = OrderedDict("cuemem" => MEM)

    notnan(x)       = OrderedDict(x  => x->((!).(isnan).(x)))
    minmax(x, m, M) = OrderedDict(x  => x-> x .>= m .&& x .<= M)
    max(x, M)       = OrderedDict(x  => x-> x .<= M)
    min(x, m)       = OrderedDict(x  => x-> x .>= m)

    """
    This is kind of more of a spike count filter than a cell count
    """
    cellcount       = OrderedDict([:time,:unit] =>
                          x->groupby_summary_cond(x, :unit,
                                                  x->x.count.>50, # > 50 spikes
                                                  nrow=>:count))
    spikecount = cellcount

    trajectory_diversity  = OrderedDict([:unit, :period] => 
            x->groupby_summary_cond(x, :unit,
                                    x->x.percount.>=3, # >= 3 trajectories
                                    :period =>
                                    (x->length(unique(x))) =>
                                    :percount))
    """
    fieldrate

    filter by rate of a sample in a field
    """
    function fieldcount(x, pos...; kws...)
    end


    function test_filt(spikes)
        println(all(combine(groupby(spikes[cellcount(spikes),:],:unit),nrow=>:count)[:,:count] .> 50))
        x = groupby_summary_filt(spikes, :unit, x->x.count.>50, nrow=>:count)
        println(all(combine(groupby(x,:unit),nrow=>:count)[:,:count] .> 50))
    end

    """
    Currently matches N filters with matching keys

    ... this could do a lot more, like function as a swapin for the
    actual merge method for Dicts, and search for matching keys
    """
    function filtmerge(D::AbstractDict...)
        newD = Dict{keytype(D[1])}{Any}()
        K = Tuple(keys(D[1]))[1]
        newD[K] = [d[K] for d in D]
        newD
    end

    function groupby_summary_filt(df, splitby, summary_condition, combine_args...)
        groups = groupby(df, splitby, sort=true)
        summaries = combine(groups, combine_args...)
        summaries[!,"condition"] = summary_condition(summaries)
        groups = groups[ summaries.condition ]
        return combine(groups, identity)
    end

    """
             groupby_summary_cond(df::DataFrame, splitby, summary_condition,
                                  combine_args...)

    splits(splitby) a df and filters by a summary condition.

    ### Algorithm
    after it splits, we transform the splits with combine_args and apply a
    condition to the resultant combined table to filter the original dataframe.

    The original and combined dataframe are kept in register by an index
    variable.
    """
    function groupby_summary_cond(df::DataFrame, splitby, summary_condition,
                                  combine_args...)::BitVector
        columns = names(df)
        if splitby isa Vector{Symbol} || splitby isa Symbol
            columns = Symbol.(columns)
        end
        if all(in.(splitby, [columns]))
            df.condition = BitVector(zeros(size(df,1)))
            df[!,:index] = 1:size(df,1)
            groups = groupby(df, splitby, sort=true)
            summaries = combine(groups, combine_args...)
            summaries[!,"condition"] = summary_condition(summaries)
            summaries = groupby(summaries, splitby)
            for (summary,group) in zip(summaries,groups)
                @assert summary[1,splitby] == group[1,splitby]
                if summary.condition[1]
                    group[!,:condition] .= true
                else
                    group[!,:condition] .= false
                end
            end
            return sort(combine(groups, identity), :index)[!,:condition]
        else
            return BitVector(ones(size(df,1)))
        end
    end

    # Create a set of predefined filter combinations
    function get_filters()
        initial = merge(speed_lib, spikecount)
        filters = OrderedDict{Symbol,Union{OrderedDict,Nothing}}()
        filters[:all]         = initial
        filters[:task]        = merge(filters[:all], Filt.task)
        filters[:correct]     = merge(filters[:all], Filt.correct)
        filters[:error]       = merge(filters[:all], Filt.error)
        filters[:nontask]     = merge(filters[:all], Filt.nontask)
        filters[:memory]      = merge(filters[:all], Filt.mem)
        filters[:cue]         = merge(filters[:all], Filt.cue)
        filters[:cue_correct] = merge(filters[:all], Filt.cue)
        filters[:cue_error]   = merge(filters[:all], Filt.cue, Filt.error)
        filters[:mem_correct] = merge(filters[:all], Filt.mem, Filt.correct)
        filters[:mem_error]   = merge(filters[:all], Filt.mem, Filt.error)
        filters[:none]        = nothing
        filters
    end

    """
        get_filter_req

    Gets the fields required by the set of filters to operate
    """
    function get_filter_req(F::AbstractDict)
        sets = [String.(f) for f ∈ keys(F)]
        sets = [s isa Vector ? s : [s] for s in sets]
        unique(collect(Iterators.flatten(sets)))
    end

end

