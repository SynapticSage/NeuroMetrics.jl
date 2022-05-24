module filt

using DataFrames
using DataStructures
export filters
using Base: merge

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

cellcount       = OrderedDict(All() =>
                              x->groupby_summary_cond(x, :unit,
                                                      x->x.count.>50,
                                                      nrow=>:count))
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
function groupby_summary_cond(df, splitby, summary_condition, combine_args...)
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
    filters = OrderedDict(:all => merge(speed_lib, cellcount))
    filters[:task]        = merge(filters[:all], task)
    filters[:correct]     = merge(filters[:all], correct)
    filters[:error]       = merge(filters[:all], error)
    filters[:nontask]     = merge(filters[:all], nontask)
    filters[:memory]      = merge(filters[:all], mem)
    filters[:cue]         = merge(filters[:all], cue)
    filters[:cue_correct] = merge(filters[:all], cue)
    filters[:cue_error]   = merge(filters[:all], cue, error)
    filters[:mem_correct] = merge(filters[:all], mem, correct)
    filters[:mem_error]   = merge(filters[:all], mem, error)
    filters
end

end

