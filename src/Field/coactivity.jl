module coactivity

using DataFrames
using Statistics
using Infiltrator
using ProgressMeter
using Polyester

import DIutils

    function get_coactive_events(u1::SubDataFrame, u2::SubDataFrame; 
            thresh=0.020, chunk=5000)
        events = Array{Float64}([])
        diffs  = Array{Float32}([])
        e1 = [extrema(u1.time)...]
        e2 = [extrema(u2.time)...]
        if any(DIutils.in_range(e1, e2)) && any(DIutils.in_range(e2,e1))
            for t = 1:chunk:length(u1.time)
                events_n1 = u1.time[t:min(t+chunk,length(u1.time))]
                e = extrema(events_n1)
                events_n2 = u2.time[DIutils.in_range(u2.time,
                                                   [e[1]-thresh, e[2]+thresh])]
                Δ = abs.(events_n1 .- events_n2[:,:]')
                hits = findall(Δ .< thresh)
                append!(events,
                      [mean([events_n1[h[1]], events_n2[h[2]]]) for h in hits])
                append!(diffs, Δ[hits])
            end
        end
        return events, diffs
    end

    function get_coactive_units(spikes::DataFrame)
        units = groupby(spikes,:unit)
        get_coactive_units(units)
    end

    function get_coactive_units(units::GroupedDataFrame; 
            append_to_spikes::Bool=false)

        last_unit = maximum(combine(units, :unit => first).unit_first)
        u_count = length(units)
        D = Vector{DataFrame}([])
        @showprogress for (i1, u1) in zip(1:u_count, units)
            #for u2 in units[i1+1:u_count]

            for u2 in units[i1+1:u_count]

                #@info "units" first(u1.unit) first(u2.unit)
                coactive_times, diffs = get_coactive_events(u1, u2)
                if !(isempty(coactive_times)) && !(isempty(diffs))
                    last_unit = last_unit + 1
                    area = first(u1.area) * "_" * first(u2.area)
                    d = Dict(:unit  => last_unit, 
                             :time  => coactive_times,
                             :diffs => diffs,
                             :unit1 => u1.unit[1],
                             :unit2 => u2.unit[1],
                             :area  => area,
                             :coact => true
                            )
                    try
                        d = DataFrame(d)
                        #@assert all(key->isdefined(d,key), keys(d))
                    catch
                        @warn "fuck"
                        println("problem")
                        @infiltrate
                    end
                    push!(D, d)
                end
            end
        end
        out = append_to_spikes ? combine(units, identity) : DataFrame()
        @exfiltrate
        try
            #append!(out, DataFrame(d); cols=:union)
            out = vcat(out, D...; cols=:union)
        catch
            @infiltrate
        end
        return out
    end

    function units_from_crosscorr()
    end


    """
    """
    function coactivecelldict(coact::DataFrame)::Dict
        units = unique(coact.unit)
        C = Dict{Int, Tuple{Int,Int,String}}() 
        for unit in units
            unit1, unit2, area = coact[findfirst(coact.unit .== unit),
                                 [:unit1, :unit2, :area]]
            C[unit] = (unit1, unit2, area)
        end
        C
    end

    function coactivecelltable(celldict::Dict)::DataFrame
        K = collect(keys(celldict))
        Vunit = hcat([[x,y] for (x,y,z) in values(celldict)]...)'
        Varea = [z for (x,y,z) in values(celldict)]
        D =OrderedDict(:unit => K, :unit1 => Vunit[:,1], :unit2 => Vunit[:,2],
                    :area => area)
        DataFrame(D)
    end
    coactivecelltable(coact::DataFrame)::DataFrame = 
        coactivecelltable(coactivecelldict(coact))


end

