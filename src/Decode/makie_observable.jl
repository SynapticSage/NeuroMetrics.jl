module makie_observable

    using DataFramesMeta
    export select_range, select_est_range, select_time, 
           select_prob, select_prob4

    function select_range(t, T; data::DataFrame=DataFrame(), Δ_bounds=nothing)
        time  = T[t]
        data = @subset(data,   (:time .> (time - Δ_bounds[1])) .&&
                               (:time .< (time + Δ_bounds[2])))
        data.time = data.time .- T[t]
        data
    end
    function select_est_range(t, T, tr, Δt, Δi; data::DataFrame=DataFrame(), Δ_bounds=nothing)
        tt = tr + (t-1)*Δi
        Δ = -Int(round(Δ_bounds[1]/Δt)) : Int(round(Δ_bounds[2]/Δt))
        start, stop = max(1,tt+Δ[1]), min(length(T),tt+Δ[2])
        center_time = data.time[t]
        data = data[start:stop,:]
        data.time .-= center_time
        data
    end
    function select_est_range(t, T, tr, Δt; data::DataFrame=DataFrame(), Δ_bounds=nothing)
        I = utils.searchsortednearest(data.time, T[t])
        if I != 1 && !(isnan(I))
            Δ = -Int(round(Δ_bounds[1]/Δt)) : Int(round(Δ_bounds[2]/Δt))
            center_time = data.time[I]
            @debug "I=$I, Δ=$Δ"
            data = data[min.(max.(I.+Δ,1), length(T)),:]
            data.time .-= center_time
        else
            data = DataFrame(lfp[1,:])
            data.time .= NaN
        end
        return data
    end
    function select_time(t, T; data::DataFrame=DataFrame(), Δ_bounds=nothing)
        time  = T[t]
        I = utils.searchsortednearest(beh.time, time)
        data = data[I, :]
        data.time = 0
        data
    end
    function select_prob(t, T; prob::Array{<:Real,3}=nothing)
        time  = min(max(t, 1), length(T))
        D = prob[:,:, time]
    end
    function select_prob4(t, T; prob::Array{<:Real,4}=nothing)
        time  = min(max(t, 1), length(T))
        D = prob[:,:, time, :]
    end

    function select_events(t, T; events::DataFrame=DataFrame())
        select = (T[t] .>= events.start) .&& (T[t] .< events.stop)
        events[select, :]
    end

end
