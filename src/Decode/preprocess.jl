# -------------------------------------------------
# PREPROCESS LFP 
# TODO Turn this into a function
# -------------------------------------------------

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

    cycles = filter(:amp_mean      => amp->(amp .> 50) .& (amp .< 600),    cycles)
    cycles = filter(:δ             => dur->(dur .> 0.025) .& (dur .< 0.5), cycles)
    cycles = filter(:velVec_median => (𝒱  -> abs.(𝒱)  .> 2),               cycles)

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
                               :start => (t->matchdxy.(t))   => [:start_x_dec, :start_y_dec],
                               :stop  => (t->matchdxy.(t))   => [:stop_x_dec,  :stop_y_dec])
    cycles = transform(cycles, [:start_x, :stop_x] => diff => :Δx,
                               [:start_y, :stop_y] => diff => :Δy,
                               [:start_y_dec, :stop_y_dec] => diff => :Δy_dec,
                               [:start_x_dec, :stop_x_dec] => diff => :Δx_dec)

    ripples = transform(ripples, :start => (t->match.(t, "x")) => :start_x,
                                 :stop  => (t->match.(t, "x")) => :stop_x,
                                 :start => (t->match.(t, "y")) => :start_y,
                                 :stop  => (t->match.(t, "y")) => :stop_y,
                                 :start => (t->matchdxy.(t))   => [:start_x_dec, :start_y_dec],
                                 :stop  => (t->matchdxy.(t))   => [:stop_x_dec,  :stop_y_dec])
    ripples = transform(ripples, [:start_x, :stop_x] => diff => :Δx,
                                 [:start_y, :stop_y] => diff => :Δy,
                                 [:start_y_dec, :stop_y_dec] => diff => :Δy_dec,
                                 [:start_x_dec, :stop_x_dec] => diff => :Δx_dec)
    ripples, cycles
end
