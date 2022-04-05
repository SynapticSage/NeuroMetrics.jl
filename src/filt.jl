module filt

using DataFrames

speed = Dict("velVec"=>x->(abs.(x) .> 4))
speed_lib = Dict("velVec"=>x->(abs.(x) .> 2))
still = Dict("velVec"=>x->(abs.(x) .< 0.5))
correct = Dict("correct" => x-> x.==1)

notnan(x) = Dict(x  => x->((!).(isnan).(x)))
minmax(x, m, M) = Dict(x  => x-> x .>= m .&& x .<= M)
max(x, M) = Dict(x  => x-> x .<= M)
min(x, m) = Dict(x  => x-> x .>= m)

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

cellcount = Dict( All() => x->groupby_summary_cond(x, :unit, x->x.count.>50, nrow=>:count))

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
function merge(D::Dict...)
    newD = Dict{keytype(D[1])}{Any}()
    K = Tuple(keys(D[1]))[1]
    newD[K] = [d[K] for d in D]
    newD
end

end

