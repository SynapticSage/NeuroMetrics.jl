export velocity_filter_ripples, get_theta_cycles, 
       curate_lfp_theta_cycle_and_phase, annotate_ripples_to_lfp,
       annotate_vector_info


function velocity_filter_ripples(beh, ripples)
    beh, lfp = raw.register(beh, lfp; transfer=["velVec"], on="time")
    beh, ripples = raw.register(beh, ripples; transfer=["velVec"], on="time")
    ripples = ripples[abs.(ripples.velVec) .< 2, :]
end

function get_theta_cycles(lfp)
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

function curate_lfp_theta_cycle_and_phase(lfp, cycles)
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
function annotate_ripples_to_lfp(lfp, ripples)
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
function annotate_vector_info(ripples, cycles)

    # Match functions
    match(time, col) = beh[utils.searchsortednearest(beh.time, time),col]
    function matchdxy(time::Real) 
        #@infiltrate
        I =  utils.searchsortednearest(T, time)
        D = replace(dat[:,:,I], NaN=>0)
        xi = argmax(maximum(D, dims=2), dims=1)
        yi = argmax(utils.squeeze(maximum(D, dims=1)), dims=1)
        [x[xi][1], y[yi][1]]
    end

    cycles = transform(cycles, :start => (t->match.(t, "x")) => :start_x,
                               :stop  => (t->match.(t, "x")) => :stop_x,
                               :start => (t->match.(t, "y")) => :start_y,
                               :stop  => (t->match.(t, "y")) => :stop_y,
                               :start => (x->matchdxy.(x))   => [:start_x_dec, :start_y_dec],
                               :stop  => (x->matchdxy.(x))   => [:stop_x_dec, :stop_y_dec])

    ripples = transform(ripples, :start => (t->match.(t, "x")) => :start_x,
                                 :stop  => (t->match.(t, "x")) => :stop_x,
                                 :start => (t->match.(t, "y")) => :start_y,
                                 :stop  => (t->match.(t, "y")) => :stop_y,
                                 :start => (x->matchdxy.(x))   => [:start_x_dec, :start_y_dec],
                                 :stop  => (x->matchdxy.(x))   => [:stop_x_dec, :stop_y_dec])
    ripples, cycles
end


function separate_theta_ripple_and_non_decodes(lfp, dat; doRipplePhase::bool=false)
    lfp = sort(combine(lfp, identity), :time)
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

function convert_to_sweeps(lfp, theta, ripple; doRipplePhase::bool=false)
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
