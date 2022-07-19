module coactivity

using DataFrames
using Statistics
using Infiltrator
using ProgressMeter
using Polyester

import Utils

    function get_coactive_events(u1::SubDataFrame, u2::SubDataFrame; 
            thresh=0.020, chunk=5000)
        events = Array{Array{Float64}}([])
        e1 = [extrema(u1.time)...]
        e2 = [extrema(u2.time)...]
        if any(Utils.in_range(e1, e2)) && any(Utils.in_range(e2,e1))
            for t = 1:chunk:length(u1.time)
                uu1 = u1.time[t:min(t+chunk,length(u1.time))]
                e = extrema(uu1)
                uu2 = u2.time[Utils.in_range(u2.time, [e[1]-thresh, e[2]+thresh])]
                Δ = abs.(uu1 .- uu2[:,:]')
                Δ = findall(Δ .< thresh)
                push!(events,
                      [mean([uu1[δ[1]], uu2[δ[2]]]) for δ in Δ]
                     )
            end
            events = sort(vcat(events...))
        end
        return events
    end

    function make_coactive_units(spikes::DataFrame)
        units = groupby(spikes,:unit)
        make_coactive_units(units)
    end

    function make_coactive_units(units::GroupedDataFrame; 
            append_to_spikes::Bool=false)
        last_unit = maximum(combine(units, :unit => first).unit_first)
        u_count = length(units)
        D = []
        @showprogress for (i1, u1) in zip(1:u_count, units)
            Threads.@threads for u2 in units[i1+1:u_count]
                @info "units" first(u1.unit) first(u2.unit)
                coactive_times = get_coactive_events(u1, u2)
                if !(isempty(coactive_times))
                    last_unit = last_unit + 1
                    area = first(u1.area) * "_" * first(u2.area)
                    d = Dict(:unit => last_unit, 
                             :time => coactive_times,
                             :area => area
                            )
                    push!(D, d)
                end
            end
        end
        @exfiltrate
        spikes = combine(units, identity)
        spikes = DataFrame()
        if append_to_spikes
            [append!(spikes, DataFrame(d); cols=:union) for d in D]
        else
            [append!(spikes, DataFrame(d); cols=:union) for d in D]
        end
        return spikes
    end



end

