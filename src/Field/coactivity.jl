module coactivity

using DataFrames
using Statistics
using Infiltrator
using ProgressMeter

    function get_coactive_events(u1::SubDataFrame, u2::SubDataFrame; 
            thresh=0.020)
        Δ = abs.(u1.time[:,:] .- u2.time[:,:]')
        Δ = findall(Δ .< thresh)
        [mean([u1.time[δ[1]], u2.time[δ[2]]]) for δ in Δ]
    end

    function make_coactive_units(spikes::DataFrame)
        units = groupby(spikes,:unit)
        make_coactive_units(units)
    end

    function make_coactive_units(units::GroupedDataFrame)
        last_unit = maximum(combine(units, :unit => first).unit_first)
        u_count = length(units)
        D = []
        #Base.length(::DataFrame) = 1
        @infiltrate
        @showprogress for (i1, u1) in zip(1:u_count, units)
            @showprogress for u2 in units[i1:u_count]
                coactive_times = get_coactive_events(u1, u2)
                last_unit = last_unit + 1
                area = first(u1.area) * "_" * first(u2.area)
                d = Dict(:unit => last_unit, 
                         :time => coactive_times,
                         :area => area
                        )
                @infiltrate
                append!(D, d)
            end
        end
        spikes = combine(units, identity)
        [append!(spikes, DataFrame(d); cols=:union) for d in D]
        return spikes
    end



end

