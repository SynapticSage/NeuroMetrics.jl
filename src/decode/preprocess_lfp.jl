export velocity_filter_ripples, get_theta_cycles, 
       curate_lfp_theta_cycle_and_phase, annotate_ripples_to_lfp,
       annotate_vector_info
export separate_theta_ripple_and_non_decodes

using DataFrames
using LoopVectorization
using Infiltrator


function velocity_filter_ripples(beh::DataFrame, ripples::DataFrame)
    beh, lfp = raw.register(beh, lfp; transfer=["velVec"], on="time")
    beh, ripples = raw.register(beh, ripples; transfer=["velVec"], on="time")
    ripples = ripples[abs.(ripples.velVec) .< 2, :]
end

function get_theta_cycles(lfp::DataFrame)
    lfp = raw.lfp.annotate_cycles(lfp, method="peak-to-peak")
    lfp.phase = raw.lfp.phase_to_radians(lfp.phase)
    cycles = raw.lfp.get_cycle_table(lfp,
                                     :velVec => (x->median(abs.(x))) => :velVec_median;
                                     end_period=:stop)
    transform!(cycles, [:start,:end] => ((x,y)->mean([x,y])) => :time)
    cycles = filter(:amp_mean => amp->(amp .> 50) .& (amp .< 600), cycles)
    cycles = filter(:δ => dur->(dur .> 0.025) .& (dur .< 0.5), cycles) 
    cycles = filter(:velVec_median => (𝒱  -> abs.(𝒱)  .> 2) , cycles)
    lfp, cycles
end

function curate_lfp_theta_cycle_and_phase(lfp::DataFrame, cycles::DataFrame)
    # TODO Remove any cycles in side a ripple
    # Remove cycle labels in lfp of bad Θ cycles
    #lfpcyc, cyc = lfp.cycle, cycles.cycle
    lfp.chunks = Int16.(round.(((1:size(lfp,1))./1000000)))
    lfp = groupby(lfp,:chunks)
    @time Threads.@threads for group in lfp
        LoopVectorization.@avxt goodcycles = in.(group.cycle, [cycles.cycle])
        group[(!).(goodcycles), :cycle] .= -1
    end
    lfp = combine(lfp, identity)
end

# (4) Annotate ripples into lfp
function annotate_ripples_to_lfp(lfp::DataFrame, ripples::DataFrame)
    ripples.type = ripples.area .* " ripple"
    ripples.rip_id = 1:size(ripples,1)
    lfp = begin
        if "area" in names(lfp)
            lfp[!,Not("area")]
        else
            lfp
        end
    end
    lfp = raw.registerEventsToContinuous(ripples, lfp, 
                             on="time", 
                             eventStart="start", 
                             eventStop="stop", 
                             ifNonMissingAppend=true,
                             targetEltype=Dict("area"=>String),
                             transfer=["rip_id","type","area"])
    lfp.rip_phase = Float32.(combine(groupby(lfp, :rip_id, sort=false),
                                     x->1/nrow(x)*ones(nrow(x))).x1)
    lfp, ripples
end

# (5) Annotate cycles with, decode vector, next/previous goal data
function annotate_vector_info(ripples::DataFrame, cycles::DataFrame, beh::DataFrame, 
                              dat::AbstractArray, x::Vector, y::Vector, T::Vector
                              current_phase::Tuple = (-pi, 0)
                              final_phase::Tuple   = (pi-pi/10, pi))

    # Cycle time based decode values
    #vecOfRow(df) = [d[1] for d in eachcol(DataFrame(df))]
    match(time, col::Vector{Symbol}) = begin
        res = beh[utils.searchsortednearest.([beh.time], time),col]
        res = [x for x in eachrow(Matrix(res))]
    end
    function matchdxy(time::Real) 
        #@infiltrate
        I =  utils.searchsortednearest(T, time)
        D = replace(dat[:,:,I], NaN=>0)
        xi = argmax(maximum(D, dims=2), dims=1)
        yi = argmax(utils.squeeze(maximum(D, dims=1)), dims=1)
        Float32.([x[xi][1], y[yi][1]])
    end
    # Phase based decode values
    function get_phase_range_start_stop(event, lfp, ϕ₀, ϕ₁)
        inds = lfp.time .>= event.start .&& lfp.time .< event.stop
        lfp = lfp[inds,:]
        start, stop = findfirst(lfp.phase .>= ϕ₀),
                      findfirst(lfp.phase .< ϕ₁)
        start, stop = lfp[start, :time],
                      lfp[stop, :time]
        start, stop
    end
    X, Y = ndgrid(x, y)
    function meandxy(start::Real, stop::Real)
        I₁,I₂ = utils.searchsortednearest(T, start),
                utils.searchsortednearest(T, stop)
        D = replace(dat[:,:,I₁:I₂], NaN=>0)
        sD = mean(D)
        x  = mean(X.*D)/sD
        y  = mean(Y.*D)/sD
        Float32.([x, y])
    end
    function get_mean_prss(event, lfp, ϕ₁, ϕ₂)
        meandxy(get_phase_range_start_stop(event,lfp,ϕ₁,ϕ₂)...)
    end

    removal_list = [:act₀, :act₁, :dec₀, :dec₁, :act₀₁, :dec₀₁]
    remove = [x for x in removal_list if x in propertynames(cycles)]
    @debug remove
    cylces = cycles[!, Not(remove)]
    cycles = transform(cycles, :start => (t->(match(t, [:x,:y]))) => :act₀,
                               :stop  => (t->(match(t, [:x,:y]))) => :act₁,
                               :start => (x->(matchdxy.(x)))   => :dec₀,
                               :stop  => (x->(matchdxy.(x)))   => :dec₁)

    cycles = transform(cycles, [:act₀,:act₁]   => ((a,b) -> b .- a) => :act₀₁,
                               [:dec₀, :dec₁]  => ((a,b) -> b .- a) => :dec₀₁)

    lfp = groupby(lfp,:cycle)
    cycles[!,:dec_ϕu]   = [Float32.([NaN, NaN]) for i in 1:size(cycles,1)]
    cycles[!,:dec_ϕd]   = cycles[:, :dec_ϕu]
    for (lf, cycle) in zip(lfp, eachrow(cycles))
        lower, upper = extrema(lf.phase)
        if !(isapprox(lower,-pi, atol=0.8)) ||
           !(isapprox(upper, pi, atol=0.8))
           continue
       end
       cycle.dec_ϕd = get_mean_prss(cycle, lf, current_phase...)
       cycle.dec_ϕu  = get_mean_prss(cycle, lf, final_phase...)
    end
    cycles[!,:dec_ϕdu] = cycles.dec_ϕu - cycles.dec_ϕd
    #cycles[!,:curfinal] = cycles.curfinal_Δx + (cycles.curfinal_Δy)im


    removal_list = [:act₀, :act₁, :dec₀, :dec₁, :act₀₁, :dec₀₁]
    remove = [x for x in removal_list if x in propertynames(ripples)]
    @debug remove
    cylces = ripples[!, Not(remove)]
    ripples = transform(ripples, :start => (t->(match(t, [:x,:y]))) => :act₀,
                               :stop  => (t->(match(t, [:x,:y]))) => :act₁,
                               :start => (x->(matchdxy.(x)))   => :dec₀,
                               :stop  => (x->(matchdxy.(x)))   => :dec₁)

    ripples = transform(ripples, [:act₀,:act₁]   => ((a,b) -> b .- a) => :act₀₁,
                               [:dec₀, :dec₁]  => ((a,b) -> b .- a) => :dec₀₁)


    ripples, cycles
end


function separate_theta_ripple_and_non_decodes(T, lfp, dat; doRipplePhase::Bool=false)
    if lfp isa GroupedDataFrame
        lfp = combine(lfp, identity)
    end
    lfp = sort(lfp, :time)
    # Theta : Create probability chunks by phase
    dat = Float32.(dat)
    theta, ripple, non = copy(dat), repeat(copy(dat), outer=(1,1,1,4)), copy(dat)
    @time Threads.@threads for (t,time) in collect(enumerate(T))
        I = utils.searchsortednearest(lfp.time, time)
        θ, ρ, n  = view(theta, :, :, t), view(ripple, :, :, t, :),
                   view(non, :, :, t)
        not_a_theta_cycle = lfp.cycle[I] == -1
        is_a_ripple = !(ismissing(lfp.rip_id[I]))
        if not_a_theta_cycle # NOT THETA
            θ .= NaN
            if is_a_ripple # IS RIPPLE?
                n .= NaN
                if lfp.area[I] == "CA1"
                    ρ[:, :,  [2, 3, 4]] .= NaN
                elseif lfp.area[I] == "CA1PFC"
                    ρ[:, :,  [1, 3, 4]] .= NaN
                elseif lfp.area[I] == "PFCCA1"
                    ρ[:, :,  [1, 2, 4]] .= NaN
                elseif lfp.area[I] == "PFC"
                    ρ[:, :,  [1, 2, 3]] .= NaN
                else
                    @error "Not a valid type"
                end

                if doRipplePhase
                    ρ[(!).(isnan.(ρ))]   .= lfp.rip_phase[I]
                end
            else
                ρ .= NaN
            end
        else # THETA CYCLE
            θ[(!).(isnan.(θ))] .= lfp.phase[I]
            ρ .= NaN
            n .= NaN
        end
    end
    return theta, ripple, non
end

function convert_to_sweeps(lfp::DataFrame, theta::Array, ripple::Array; 
                           doRipplePhase::Bool=false)

    # Create cumulative theta sweeps
    sweep = (a,b)->isnan(b) ? a : nanmean(cat(a, b, dims=3), dims=3)
    lfp = groupby(lfp,:cycle)
    @time Threads.@threads for group in lfp
        cycStart, cycStop = utils.searchsortednext(T, group.time[1]),
                            utils.searchsortednext(T, group.time[end])
        cycle = cycStart:cycStop
        if cycle == 1:1 || group.cycle[1] == -1 || ismissing(group.cycle[1])
            continue
        end
        θ = view(theta, :, :, cycle)
        θ = accumulate(sweep, theta[:,:,cycle],  dims=3)
    end
    lfp = combine(lfp, identity)
    lfp = sort!(lfp,:time)

    # Cumulative ripple sweeps
    if doRipplePhase
        lfp = groupby(lfp, :rip_id)
        @time @Threads.threads for group in lfp
            cycStart, cycStop = utils.searchsortednext(T, group.time[1]),
                                utils.searchsortednext(T, group.time[end])
            local cycle = cycStart:cycStop
            if cycle == 1:1 || group.cycle[1] == -1 ||
               ismissing(group.cycle[1])
                continue
            end
            ripple[:, :, cycle] = accumulate((a,b)->sweep.(a,b),
                                             ripple[:,:,cycle], dims=3)
        end
        lfp = combine(lfp,identity)
    end
    return theta, ripple

end

function beh_to_cycles(beh, cycles, cycreg=:time, behreg=:time;
        register=[:stopWell,:futureStopWell,:pastStopWell,:relativetraj])
    @assert cycreg == behreg "For now these have to be the same...."
end

function annotate_behavior_to_cycles(beh::DataFrame, 
        events::DataFrame, pertrajlabel=:traj)
    if :time ∉ propertynames(events)
        events[!,:time] = vec(mean([events.start events.end],dims=2))
    end
    if :cycle ∈ propertynames(events)
        cycle_unit = :cycle
    elseif :rip_id ∈ propertynames(events)
        cycle_unit = :rip_id
    end
    transfer = String.([:traj, :correct, :stopWell, :futureStopWell, :pastStopWell,
                        :stopWell])
    _, events = raw.register(beh, events, on="time", transfer=transfer)
    groups = groupby(events, :traj)
    for group in groups
        group[!,:cycle_traj] = replace(group[!,cycle_unit],-1=>missing)
        group[!,:cycle_traj] = group[!,:cycle_traj] .- minimum(group[!,:cycle_traj]) .+ 1
    end
    events = combine(groups, identity)
    x = events[!,:cycle_traj]
    events[!,:cycle_traj] = convert(Vector{Float32}, coalesce(x, missing=>NaN))
end

function annotate_explodable_cycle_metrics(beh::DataFrame, 
        events::DataFrame, dat::AbstractArray,
        x::Vector{<:Real}, y::Vector{<:Real}, T::Vector{<:Real}
    )

    # Cycle time based decode values
    function imatchdxy(I::Real) 
        #@infiltrate
        D = replace(dat[:,:,I], NaN=>0)
        xi = argmax(maximum(D, dims=2), dims=1)
        yi = argmax(utils.squeeze(maximum(D, dims=1)), dims=1)
        Float32.([x[xi][1], y[yi][1]])
    end

    events.midpoint = vec(mean([events.start events.stop],dims=2))
    events.time     = events.midpoint
    _, events = raw.register(beh,events,on="time",transfer=["traj"])

    #c⃗ ᵢⱼ, trajreltime, time
    events.cij         = Vector{Vector}(undef,size(events,1))
    events.trajreltime = Vector{Vector}(undef,size(events,1))
    events.time        = Vector{Vector}(undef,size(events,1))
    events.trajtime    = Vector{Vector}(undef,size(events,1))
    P = Progress(size(events,1), dt=0.1, 
                 desc="Adding explodable fields to events")
    Threads.@threads for row in eachrow(events)
        if row.cycle == -1
            row.cij         = []
            row.time        = []
            row.trajreltime = []
        end
        Tind = findall(T .>= row.start .&& T .< row.stop)
        row.cij = imatchdxy.(Tind)
        #row.cij_x, row.cij_y = tmp1, tmp2
        row.time = T[Tind]
        row.trajtime = T[Tind] .- minimum(T[Tind])
        row.trajreltime = utils.searchsortednearest.([beh.time], T[Tind])
        row.trajreltime = beh[row.trajreltime, :trajreltime]
        cumchange = mean(cumsum(diff(row.trajreltime)))
        if cumchange > 0
            row.trajreltime = row.trajreltime[begin]:cumchange:row.trajreltime[end]
        end
        next!(P)
    end
    # trajcycletime
    events = events[events.cycle .!=-1,:]
    events = groupby(events, :traj)
    for event in events
        global prevEnd
        event.trajcycletime = event.trajtime
        prevEnd = NaN
        for (i,item) in enumerate(event.trajcycletime)
            global prevEnd
            if i > 1
                item .-= item[1] - prevEnd
            else
                item .-= item[1]
            end
            prevEnd = item[end]
        end
    end
    events = combine(events, identity)

end

explode_cols = [:cij, :trajreltime, :time, :trajtime, :trajcycletime]


"""
Remaps a pair of dataframe columns (the vector) to the coordinates of
an upcoming goal
"""
function annotate_vector_relative_to_goal()
    # Get best goal
    # Get angle relative to F₁, F₂, P₁
end
