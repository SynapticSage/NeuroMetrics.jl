export velocity_filter_ripples, get_theta_cycles, 
       curate_lfp_theta_cycle_and_phase, annotate_ripples_to_lfp,
       annotate_vector_info
export separate_theta_ripple_and_non_decodes
export velocity_filter_ripples

using DataFrames, DataFramesMeta
using LoopVectorization
using Infiltrator
using ProgressMeter
using LazyGrids: ndgrid
import .raw
import .table

to_complex(x) = ComplexF64(x...)
explode_cols = [:decᵢ, :decᵢᵢ, :trajreltime, :time, :trajtime, :trajcycletime]
reference_point = Dict(
    :dec₀₁=>:dec₀,
    :dec₀ᵢ=>:dec₀,
    :decᵢᵢ=>:decᵢ)
reference_time = Dict(
    :dec₀₁=>:start,
    :dec₀ᵢ=>:start,
    :decᵢᵢ=>:time)

function velocity_filter_ripples(beh::DataFrame, ripples::DataFrame)
    beh, ripples = raw.register(beh, ripples; transfer=["velVec"], on="time")
    ripples = ripples[abs.(ripples.velVec) .< 2, :]
end

function get_theta_cycles(lfp::DataFrame, beh::DataFrame)
    beh, lfp = raw.register(beh, lfp; transfer=["velVec"], on="time")
    lfp = raw.lfp.annotate_cycles(lfp, method="peak-to-peak")
    lfp.phase = raw.lfp.phase_to_radians(lfp.phase)
    cycles = table.get_periods(lfp, "cycle", 
                   :amp=>mean,
                   :velVec => (x->median(abs.(x))) => :velVec_median; 
                   end_period=:stop)
    transform!(cycles, [:start,:stop] => ((x,y)->mean([x,y])) => :time)
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
function annotate_vector_info(ripples::DataFrame, cycles::DataFrame,
        beh::DataFrame, lfp::DataFrame, dat::AbstractArray, x::Vector,
        y::Vector, T::Vector, current_phase::Tuple = (-pi, 0),
        final_phase::Tuple   = (pi-pi/10, pi))

    # Cycle time based decode values
    #vecOfRow(df) = [d[1] for d in eachcol(DataFrame(df))]
    match(time, col::Vector{Symbol}) = begin
        res = beh[utils.searchsortednearest.([beh.time], time),col]
        res = [ComplexF64(x[1]+x[2]im) for x in eachrow(Matrix(res))]
    end
    function matchdxy(time::Real) 
        I =  utils.searchsortednearest(T, time)
        D = replace(dat[:,:,I], NaN=>0)
        xi = argmax(maximum(D, dims=2), dims=1)
        yi = argmax(utils.squeeze(maximum(D, dims=1)), dims=1)
        ComplexF64(x[xi][1] + y[yi][1]im)
    end
    # Phase based decode values
    function get_phase_range_start_stop(event, lfp, ϕ₀, ϕ₁)
        inds = lfp.time .>= event.start .&& lfp.time .< event.stop
        lfp  = lfp[inds,:]
        start, stop = findfirst(lfp.phase .>= ϕ₀),
                      findfirst(lfp.phase .< ϕ₁)
        if start == nothing || stop == nothing
            start, stop = nothing, nothing
        else
            start, stop = lfp[start, :time],
                          lfp[stop, :time]
        end
    end
    X, Y = ndgrid(x, y)
    function meandxy(start::Real, stop::Real)
        I₁,I₂ = utils.searchsortednearest(T, start),
                utils.searchsortednearest(T, stop)
        D = replace(dat[:,:,I₁:I₂], NaN=>0)
        sD = mean(D)
        xm  = mean(X.*D)/sD
        ym  = mean(Y.*D)/sD
        ComplexF64(xm .+ (ym)im)
    end
    function get_mean_prss(event, lfp, ϕ₁, ϕ₂)
        start, stop = get_phase_range_start_stop(event,lfp,ϕ₁,ϕ₂)
        if (start,stop) == (nothing, nothing)
            ComplexF64(NaN + (NaN)im)
        else
            meandxy(start, stop)
        end
    end

    removal_list = [:act₀, :act₁, :dec₀, :dec₁, :act₀₁, :dec₀₁]
    remove       = [elem for elem in removal_list if elem in propertynames(cycles)]
    @debug remove
    cylces = cycles[!, Not(remove)]
    cycles = transform(cycles, :start => (t->(match(t, [:x,:y]))) => :act₀,
                               :stop  => (t->(match(t, [:x,:y]))) => :act₁,
                               :start => (x->(matchdxy.(x)))   => :dec₀,
                               :stop  => (x->(matchdxy.(x)))   => :dec₁)

    cycles = transform(cycles, [:act₀,:act₁]   => ((a,b) -> b .- a) => :act₀₁,
                               [:dec₀, :dec₁]  => ((a,b) -> b .- a) => :dec₀₁)

    lfp = groupby(lfp,:cycle)
    cycles[!,:dec_ϕu]   = [ComplexF64(NaN + (NaN)im) for i in 1:size(cycles,1)]
    cycles[!,:dec_ϕd]   = cycles[:, :dec_ϕu]
    for (lf, cycle) in zip(lfp, eachrow(cycles))
        lower, upper = extrema(lf.phase)
        if !(isapprox(lower,-pi, atol=0.8)) ||
           !(isapprox(upper, pi, atol=0.8))
           continue
       end
       cycle.dec_ϕd  = get_mean_prss(cycle, lf, current_phase...)
       cycle.dec_ϕu  = get_mean_prss(cycle, lf, final_phase...)
    end
    cycles[!,:dec_ϕdu] = cycles.dec_ϕu - cycles.dec_ϕd


    removal_list = [:act₀, :act₁, :dec₀, :dec₁, :act₀₁, :dec₀₁]
    remove = [elem for elem in removal_list if elem in propertynames(ripples)]
    @debug remove
    cylces = ripples[!, Not(remove)]
    ripples = transform(ripples, :start => (t->(match(t, [:x,:y]))) => :act₀,
                                 :stop  => (t->(match(t, [:x,:y]))) => :act₁)
    ripples = transform(ripples, 
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
        E::DataFrame, pertrajlabel=:traj)
    if :time ∉ propertynames(E)
        E[!,:time] = vec(mean([E.start E.end],dims=2))
    end
    if :cycle ∈ propertynames(E)
        cycle_unit = :cycle
    elseif :rip_id ∈ propertynames(E)
        cycle_unit = :rip_id
    end
    transfer = String.([:traj, :correct, :stopWell, :futureStopWell, :pastStopWell,
                        :stopWell])
    _, E = raw.register(beh, E, on="time", transfer=transfer)
    groups = groupby(E, :traj)
    for group in groups
        group[!,:cycle_traj] = replace(group[!,cycle_unit],-1=>missing)
        group[!,:cycle_traj] = group[!,:cycle_traj] .- minimum(group[!,:cycle_traj]) .+ 1
    end
    E = combine(groups, identity)
    x = E[!,:cycle_traj]
    E[!,:cycle_traj] = convert(Vector{Float32}, coalesce(x, missing=>NaN))
end

function annotate_explodable_cycle_metrics(beh::DataFrame, 
        E::DataFrame, dat::AbstractArray,
        x::Vector{<:Real}, y::Vector{<:Real}, T::Vector{<:Real}
    )

    if :trajreltime ∉ propertynames(beh)
        raw.behavior.annotate_relative_xtime!(beh)
    end

    if :cycle in propertynames(E)
        cycle_field = :cycle
    else
        cycle_field = :rip_id
    end


    # Cycle time based decode values
    function imatchdxy(I::Real) 
        D = replace(dat[:,:,I], NaN=>0)
        xi = argmax(maximum(D, dims=2), dims=1)
        yi = argmax(utils.squeeze(maximum(D, dims=1)), dims=1)
        Float32.([x[xi][1], y[yi][1]])
    end

    E.midpoint = vec(mean([E.start E.stop],dims=2))
    E.time     = E.midpoint
    _, E = raw.register(beh,E,on="time",transfer=["traj"])

    #c⃗ ᵢⱼ, trajreltime, time
    E.dec_0i        = Vector{Vector}(undef,size(E,1))
    E.dec_ii       = Vector{Vector}(undef,size(E,1))
    E.trajreltime_i = Vector{Vector}(undef,size(E,1))
    E.time_i        = Vector{Vector}(undef,size(E,1))
    E.trajtime_i    = Vector{Vector}(undef,size(E,1))
    P = Progress(size(E,1), dt=0.1, 
                 desc="Adding explodable fields to E")
    Threads.@threads for row in eachrow(E)
        if row[cycle_field] == -1
            row.decᵢ        = []
            row.decᵢᵢ       = []
            row.time        = []
            row.trajreltime = []
        end
        Tind = findall(T .>= row.start .&& T .< row.stop)
        row.dec_i  = imatchdxy.(Tind) # looks up vector at each time
        row.dec_ii = row.decᵢ[2:end] .- row.dec_i[1:end-1] # looks up vector at each time
        row.dec_ii = cat(NaN + (NaN)im, row.dec_ii; dims=1)
        row.time_i = T[Tind]
        row.trajtime_i = T[Tind] .- minimum(T[Tind])
        row.trajreltime_i = utils.searchsortednearest.([beh.time], T[Tind])
        row.trajreltime_i = beh[row.trajreltime, :trajreltime]
        cumchange = mean(cumsum(diff(row.trajreltime)))
        if cumchange > 0
            row.trajreltime = LinRange(row.trajreltime[begin],
                                       row.trajreltime[end],
                                       length(row.trajreltime))
            #TODO this may not be 100% accurate for sequences
            # that have a situation as follows. Suppose
            # [1,1,1,1,2,2,2,2,3,3,3] where when you translate these into
            # their respective fraciontions
            # [1,1.25,1.5,...,3,3.25,3.5] the top number is not equal to
            # the [end] number of the sequence, as above.
            # Hence there's a tiny fudge factor in its current form
        end
        next!(P)
    end

    # trajcycletime
    E = E[E[!,cycle_field] .!=-1,:]
    E = groupby(E, :traj)
    for event in E
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
    E = combine(E, identity)

end


"""
Remaps a pair of dataframe columns (the vector) to the coordinates of
an upcoming goal
"""
# Get best goal
# Get angle relative to F₁, F₂, P₁
function annotate_vector_relative_to_goal(beh::DataFrame, E::DataFrame,
    wells::DataFrame, vector=:dec₀₁, to_which=:stopWell)

    if to_which isa Real
        to_which = to_which * ones(size(beh,1))
    else
        to_which = beh[!,to_which]
    end

    function get_well_record(i)
        cols = [:x, :y]
        if i > 0
            w=wells[i, cols]
        else
            w = wells[1,cols]
            w .= NaN
        end
        w=DataFrame(w)
        #w[!,to_which] = [i]
        w
    end
    wellrecord      = vcat(get_well_record.(to_which)...)
    wellrecord.time = beh.time

    # Do we need to explode fields?
    vecofvec = typeof.(eachcol(E[!,explode_cols])) .<: 
        Vector{T} where T <: Union{Vector, Union{Missing, Vector}}
    if vector ∈ explode_cols && any(vecofvec)
        E = flatten(E, explode_cols[vecofvec])
    end


    # Start of the decode change vector, the reference point
    reference_vectors = cat(E[!, reference_point[vector]]...,dims=2)'
    # Well record registered to the decode time
    registration = utils.searchsortednearest.([wellrecord.time],
                                             E[!,reference_time[vector]]);
    wellrec =  Matrix(wellrecord[registration,[:x,:y]])
    # Vector from the well to the reference point
    point_to_well_vectors = wellrec - reference_vectors

    # Decode vectors
    decode_vectors = cat(E[!, vector]...,dims=2)'
    unitvec(x) = x./abs(x)
    decode_to_well = unitvec(decode_vectors).-unitvec(point_to_well_vectors)

    # And now we measure consistency of our decode vector to it
    abs.(decode_vectors) * angle.(decode_to_well)
end

function clean_vec_fields(events::DataFrame)
    fields = [field for field in names(events)
              if occursin("_x",field) || occursin("_y",field)]
    events = events[!, Not(fields)]
end

"""

Params
------
handle_multiple_tets = :and  | :or
:and sums correlation of all pairs, and uses that as a metric
:or sums correlation of all pairs and takes the highest correlation
"""
function add_correlation_coordination(L::DataFrame;
        area1=[], area2=[], shifts=[0], deltaB=0.3, 
		metrics=[:summary, :each])::DataFrame

	L = sort(L[:, [:time,:tetrode,:raw,:cycle]], [:time, :tetrode])
	inds = utils.squeeze(any(Int64.(L.tetrode) .∈ [area1... area2...], dims=2))
    L = L[inds, :]
    ulfp = unstack(L, :time, :tetrode, :raw)
	ulfp.cycle = @subset(lfp, :tetrode .== 5).cycle

	dT = median(diff(ulfp.time))
	w = Int(round(deltaB/dT))

	M = Matrix(ulfp[!,2:end])
	M = convert.(eltype(M), M)
	area1_locs = 1:length(area1)
	area2_locs = maximum(area1_locs)+1:maximum(area1_locs)+length(area2)
	shifts = -100:5:100 # tidbit more than 0.05 second before and after
	shifts = -200:5:200 # tidbit more than 0.1 second before and after
	shifts_num = shifts .* (median(diff(ulfp.time)))

	C = zeros(length(shifts), size(M,1), length(area1_locs) * length(area2_locs))
	@time Threads.@threads for i in 1:size(ulfp,1)
		for (s,shift) in enumerate(shifts)
			if shift == 0
				 # We're removing zero phase lag effects here
				continue
			end
			iL1, iH1 = max(i-w, 1), min(i+w, size(ulfp,1))
			iL2, iH2 = max(i-w+shift, 1), min(i+w+shift, size(ulfp,1))
			M1 = view(M, iL1:iH1,:)
			M2 = view(M, iL2:iH2,:)

			m1, m2 = M1[:, area1_locs], M2[:, area2_locs]
			if size(m1,1) != size(m2,1)
				continue
			end
			c  = m1' * m2
			c  = vec(c)
			cc = view(C,s,i,:)
			cc[:] = c
		end
	end

	# Get dynamic measurements (best shift at each time)
	# --------------------------------------------------
	dynamic_C_shiftOpt = utils.squeeze(maximum(C, dims=1))
	dynamic_S_shiftOpt = utils.squeeze(argmax(C, dims=1))
	dynamic_S_shiftOpt = [S_shiftOpt[i,j].I[1] 
								for i in 1:size(S_shiftOpt,1), 
									j in 1:size(S_shiftOpt,2)]
	dynamic_S_shiftOpt = shifts_num[S_shiftOpt]

	# Get static measurements (best shift at each time)
	# --------------------------------------------------
	indcycle = [1:size(ulfp,1) ulfp.cycle]
	for cyc in unique(indcycle[:,2])
		inds = indcycle[indcycle[:,2] .== cyc, 1]
		#Slice
		c = view(C, :, inds, :)
		#get best shift	
		# TODO
	end
	
	results = DataFrame()
	for metric in metrics
		if metric == :summary
		elseif metric == :each
			for col in eachcol(C)
			end
		end
	end

	return results
end
function add_correlation_coordination(L1::DataFrame, L2::DataFrame;
        kws...)::DataFrame
	L = vcat(L1, L2, cols=:intersect)	
	area1, area2 = unique(L1.tetrode), unique(L2.tetrode)
	add_correlation_coordination(L; area1=area1, area2=area2,
								 kws...)
end
