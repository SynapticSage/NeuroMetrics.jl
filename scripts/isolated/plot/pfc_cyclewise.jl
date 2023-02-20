
# READY DATATYPES AND MODULES
if !(:lfp in names(Main))

    # IMPORTS AND DATA
    # --------
    using GLM, Lasso, Distributions, ThreadSafeDicts
    @time include(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")

    # CONSTANTS and trackers
    # --------
    val = :value
    opt["matched"] = 3    # how many matched non-iso cycles to add per iso cycle
    opt["has_df"] = false # tracks state about whether df var for glm is backed up in a jld2 file
    props = [:x, :y, :speed, :startWell, :stopWell]

    # DATAFRAME-IZE firing rate matrix
    # and add properties of interest
    # --------------------------------
    # Get dataframe of R
    Rdf = DataFrame(R; name=val)
    # Set cycle via start cycle times
    cycles.time = cycles.start
    cycles.cycle = 1:size(cycles,1)
    DIutils.filtreg.register(cycles, Rdf, on="time", 
                             transfer=["cycle"], 
                             match=:prev)
    DIutils.filtreg.register(beh, Rdf, on="time", 
                            transfer=String.(props))

    # Figure out cycles with some number of isolated spikes
    begin
        cycstat = combine(groupby(@subset(spikes, :area .== "CA1"), 
                        :cycle), 
                :isolated => sum, 
                :unit => (n->unique(n)) => :diversity,
                [:isolated, :unit] => ((i,u)->(length(unique(u[i.==true])))) => :isodiversity,
                [:isolated, :unit] => ((i,u)->(length(unique(u[i.==false])))) => :adjdiversity,
               )
        dropmissing!(cycstat, :cycle)
    end

    # Register datasets
    @showprogress "registration" for (name,obj) in zip(("cycle","spikes","Rdf"),
                                                       (cycles, spikes, Rdf))
        @info "register" name
        DIutils.filtreg.register(cycstat, obj, on="cycle", 
                         transfer=String.(setdiff(propertynames(cycstat),[:cycle])))
    end
    DIutils.filtreg.register(cells, Rdf, on="unit", transfer=["area"])

    # Annotate spikes and Rdf
    dropmissing!(Rdf, :isolated_sum)
    Rdf_cycles    = groupby(Rdf, [:cycle])
    indexers = [:time,:isolated_sum]

    begin
        kws=(;label="")
        plot(
             (@df cycstat histogram(:isolated_sum ;xlabel="iso spikes emitted", kws...)),
             (@df cycstat histogram(:isodiversity;xlabel="# of uniq iso cells", kws...)),
             (@df cycstat histogram(:adjdiversity;ylabel="# of uniq adj cells", kws...)),
             Plot.blank((plot();Plots.annotate!(0.5,0.5,text("Cycle statistics",14))), visible=false, size=(100,50)),
             layout=grid(2,2))
    end
end


# ----------------------------------------
# Obtain adaptive grid and match θ cycles
# ----------------------------------------
# ( For matching isolated theta cycles)
begin

    # Obtain average of properties of interest per θ cycle
    # -----------------------------------------------------
    for prop in props
        cycles[!,prop] = Vector{Union{Missing,eltype(beh[!,prop])}}(missing, size(cycles,1))
    end
    E = Threads.Atomic{Int}(0)
    Threads.@threads for i in collect(1:size(cycles,1))
        cycle = cycles[i,:]
        inds = DIutils.in_range(beh.time, [cycle.start, cycle.stop])
        for prop in props
            try
                cycle[prop] = nanmean(collect(skipmissing(beh[inds, prop])))
            catch
                Threads.atomic_add!(E, 1)
            end
        end
    end
    println("Fail percent ", E[] / size(cycles,1))
    for prop in props
        cycles[!,prop] = replace(cycles[!,prop], NaN=>missing)
    end

    # Get the grid binning and occupancy of these cycles
    # -----------------------------------------------------
    speed_1σ  = Float32(std(beh.speed))
    theta_cycle_req = 1f0 # disable required number of cy
    theta_cycle_dt = Float32(median(diff(cycles.start)))
    use_behavior = false
    epsilon = 1f-6
    widths  = [10f0,  10f0,  speed_1σ, 1f0,   1f0]
    grid_kws =
            (;widths,
              radiusinc    = [0.2f0, 0.2f0, 0f0,      0f0,   0f0],
              maxrad       = [6f0,   6f0,   0.4f0,    0.4f0, 0.4f0],
              radiidefault = widths./2 .- epsilon,
              steplimit=2,
              dt=use_behavior ? Float32(median(diff(beh))) : theta_cycle_dt,
              thresh=theta_cycle_req
             )
    grd = DIutils.binning.get_grid(use_behavior ? beh : cycles, props; grid_kws...)
    occ = DIutils.binning.get_occupancy_indexed(cycles, grd)
    println("Percent cycles counted ", sum(occ.count)/size(dropmissing(cycles,props),1))

    # Save!!!
    # -----------------------------------------------------
    has_df = false
    jldopen(path_iso(opt; append="_cyclewise"), "a") do storage
        if "grd" in keys(storage)
            delete!(storage, "grd")
        end
        storage["grd"] = grd
        if "occ" in keys(storage)
            delete!(storage, "occ")
        end
        storage["occ"] = occ
        if "cycles" in keys(storage)
            delete!(storage, "cycles")
        end
        storage["cycles"] = cycles
        has_df = "df" in keys(storage)
    end

    # Match the cycles
    cycles.hasocc = (!).(ismissing.(occ.datainds))
    DIutils.filtreg.register(cycles, Rdf, transfer=["hasocc"], on="cycle")
    iso_cycles = unique(@subset(Rdf, :isolated_sum .> 0, :hasocc .== true).cycle)
    cycles.matched = Vector{Union{Vector{Int32}, Missing}}(missing, size(cycles,1))
    Threads.@threads for cyc in iso_cycles
        poss = [] 
        # Lookup cycles that match this isolated spike cycle's animal behavior
        for gridmatch in occ.datainds[cyc], cycmatch in occ.inds[gridmatch]
            push!(poss, cycmatch)
        end
        # Lookup which cycles lack isolated spikes
        poss = poss[cycles[poss, :].isolated_sum .=== 0]
        samples_to_grab = min(length(poss), opt["matched"])
        if samples_to_grab > 0
            cycles[cyc,:].matched = sample(poss, samples_to_grab, replace=false)
        end
    end
end

# GET CYCLEWISE INFORMATION
fn = path_iso(opt; append="_cyclewise")
if (!isfile(fn) && has_df) || opt["overwrite"]

    df = Vector{Any}(undef, 16)
    for thread in 1:Threads.nthreads()
        df[thread] = Vector{Union{Missing,DataFrame}}(missing, length(iso_cycles))
    end
    matched_cycle_holder = Vector{Union{Int,Missing}}(missing, opt["matched"])
    cyc_error =  Dict() 
    Infiltrator.clear_disabled!()
    prog = Progress(length(iso_cycles); desc="cycle df")
    Threads.@threads for (i,cyc) in collect(enumerate(iso_cycles))
        try
            tid = Threads.threadid()
            # Push the isolated cycle and its preceding following cycles
            push!(df[tid], grab_cycle_data(Rdf_cycles, cyc; val, indexers))
            matched_cycs = @subset(cycles, :cycle .== cyc).matched[1]
            # Push MATCHED cycles
            for mc in matched_cycs
                push!(df[tid], grab_cycle_data(Rdf_cycles, mc; val, indexers))
            end
            next!(prog)
        catch exception
            cyc_error[cyc] = exception
        #     if mod(i, 100) == 0
        #         @info cyc_error
        #     end
            sleep(0.1)
            next!(prog)
        end
    end
    @info cyc_error
    vcatnonmiss(df) = vcat(df[(!).(ismissing.(df))]...)
    df = vcatnonmiss.(df)
    df = vcatnonmiss(df)
    df.has_iso = df.isolated_sum .> 0
    # Spike count
    neuroncols = names(df)[tryparse.(Int, names(df)) .!== nothing]
    df[:,neuroncols] .*= median(diff(beh.time)) # TODO not INT because it's gaussian smoothed
    df[:,neuroncols] .= round.(df[:,neuroncols])
    df = transform(df, neuroncols .=> n -> convert(Vector{Int64}, n), renamecols=false)

    # Checkpoint
    rm(fn)
    jldsave(fn; df)
else
    storage = JLD2.jldopen(fn)
    df = storage["df"]
    close(storage)
end

# Clean data frame

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF SPIKE COUNTS  🔺
# ========================

# Now let's take the formula and apply them
formulae, models   = OrderedDict(), ThreadSafeDict()
formulae["ca1pfc"] = construct_predict_spikecount(df, cells, "CA1");
formulae["pfcca1"] = construct_predict_spikecount(df, cells, "PFC");
@assert !isempty(first(values(formulae)))
glmsets = []
for indep in ("ca1pfc", "pfcca1"), relcyc in -8:8, f in formulae[indep]
    push!(glmsets, (indep, relcyc, f))
end

prog = Progress(length(glmsets); desc="GLM spike counts")
Threads.@threads for (indep, relcyc, f) in glmsets
    try
        d = @subset(df, :relcycs .== relcyc)
        y, XX = modelcols(apply_schema(f, schema(f, d)), d) 
        FittedPoisson = fit_mle(Poisson, Int.(y))
        models[(;indep, relcyc)] = m =  fit(GLM.GeneralizedLinearModel, XX, y, FittedPoisson)
    catch
        models[(;indep, relcyc)] = nothing
    end
    next!(prog)
end

@info "saving spikcount"
if isdefined(Main, :storage); close(storage); end
storage = jldopen(fn, "a")
"model_spikecount" in keys(storage) ? delete!(storage, "model_spikecount") : nothing
global storage["model_spikecount"] = model_spikecount = models
@info "saved"
if isdefined(Main, :storage); close(storage); end

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF HAS_ISO
# ========================

# Now let's take the formula and apply them
formulae, models   = OrderedDict(), ThreadSafeDict()
formulae["ca1pfc"] = construct_predict_iso(df, cells, "CA1", :has);
formulae["pfcca1"] = construct_predict_iso(df, cells, "PFC", :has);
@assert !isempty(first(values(formulae)))
glmsets = []
for indep in ("ca1pfc", "pfcca1"), relcyc in -8:8, f in formulae[indep]
    push!(glmsets, (indep, relcyc, f))
end

prog = Progress(length(glmsets); desc="GLM has iso")
Threads.@threads for (indep, relcyc, f) in glmsets
    try
        d = @subset(df, :relcycs .== relcyc)
        y, X = modelcols(apply_schema(f, schema(f, d)), d) 
        FittedBinomial = fit_mle(Binomial, 1, y)
        models[(;indep, relcyc)] = glm(X, y, FittedBinomial)
    catch
        models[(;indep, relcyc)] = nothing
    end
    next!(prog)
end

@info "saving spikcount"
if isdefined(Main, :storage); close(storage); end
storage = jldopen(fn, "a")
"model_hasiso" in keys(storage) ? delete!(storage, "model_hasiso") : nothing
global storage["model_hasiso"] = model_hasiso = models
@info "saved"
if isdefined(Main, :storage); close(storage); end

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF ISO COUNT
# ========================

# Now let's take the formula and apply them
formulae, models   = OrderedDict(), ThreadSafeDict()
formulae["ca1pfc"] = construct_predict_iso(df, cells, "CA1", :count);
formulae["pfcca1"] = construct_predict_iso(df, cells, "PFC", :count);
@assert !isempty(first(values(formulae)))
glmsets = []
for indep in ("ca1pfc", "pfcca1"), relcyc in -8:8, f in formulae[indep]
    push!(glmsets, (indep, relcyc, f))
end

prog = Progress(length(glmsets); desc="GLM iso counts")
Threads.@threads for (indep, relcyc, f) in glmsets
    try
        d = @subset(df, :relcycs .== relcyc)
        y, XX = modelcols(apply_schema(f, schema(f, d)), d) 
        FittedPoisson = fit_mle(Poisson, Int.(y))
        models[(;indep, relcyc)] = m =  fit(GLM.GeneralizedLinearModel, XX, y, FittedPoisson)
    catch exception
        models[(;indep, relcyc)] = exception
    end
    next!(prog)
end

@info "saving isocount"
if isdefined(Main, :storage); close(storage); end
storage = jldopen(fn,"a")
"model_isocount" in keys(storage) ? delete!(storage, "model_isocount") : nothing
global storage["model_isocount"] = model_isocount = models
@info "saved"
if isdefined(Main, :storage); close(storage); end


#    _  _     ____        _         __                     _             
#  _| || |_  |  _ \  __ _| |_ __ _ / _|_ __ __ _ _ __ ___ (_)_ __   __ _ 
# |_  ..  _| | | | |/ _` | __/ _` | |_| '__/ _` | '_ ` _ \| | '_ \ / _` |
# |_      _| | |_| | (_| | || (_| |  _| | | (_| | | | | | | | | | | (_| |
#   |_||_|   |____/ \__,_|\__\__,_|_| |_|  \__,_|_| |_| |_|_|_| |_|\__, |
#                                                                  |___/ 
func = v->if v isa Exception
    NaN
else
    adjr2(v, :devianceratio)
end

plot(
    glmplot(model_spikecount, "pfcca1", func, label="spike count, general"),
    glmplot(model_isocount,    "pfcca1", func, label="iso count"),
    glmplot(model_hasiso,     "pfcca1", func, label="iso occurance");
    ylims=(0,1),
    link=:y,
    layout=grid(3,1)
)


plot(
    glmplot(model_spikecount, "ca1pfc", func, label="spike count, general"),
    glmplot(model_isocount,    "ca1pfc", func, label="iso count"),
    glmplot(model_hasiso,     "ca1pfc", func, label="iso occurance");
    ylims=(0,1),
    link=:y,
    layout=grid(3,1)
)
