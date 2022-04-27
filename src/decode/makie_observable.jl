
function select_range(t, data=spikes, Δ_bounds=Δ_bounds)
    time  = T[t]
    data = @subset(data,   (:time .> (time - Δ_bounds[1])) .&&
                           (:time .< (time + Δ_bounds[2])))
    data.time = data.time .- T[t]
    data
end
function select_est_range(t, tr, Δt, Δi, data=beh, Δ_bounds=Δ_bounds)
    tt = tr + (t-1)*Δi
    Δ = -Int(round(Δ_bounds[1]/Δt)) : Int(round(Δ_bounds[2]/Δt))
    start, stop = max(1,tt+Δ[1]), min(length(T),tt+Δ[2])
    center_time = data.time[t]
    data = data[start:stop,:]
    data.time .-= center_time
    data
end
function select_est_range(t, tr, Δt, data=beh, Δ_bounds=Δ_bounds)
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
function select_time(t, data=spikes, Δ_bounds=Δ_bounds)
    time  = T[t]
    I = utils.searchsortednearest(beh.time, time)
    data = data[I, :]
    data.time = 0
    data
end
function select_prob(t, prob::Array{<:Real,3}=dat)
    time  = min(max(t, 1), length(T))
    D = prob[:,:, time]
end
function select_prob4(t, prob::Array{<:Real,4}=dat)
    time  = min(max(t, 1), length(T))
    D = prob[:,:, time, :]
end
