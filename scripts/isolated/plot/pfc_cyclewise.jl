
# READY DATATYPES AND MODULES
if !(:lfp in names(Main))

    # IMPORTS AND DATA
    # --------
    using GLM, Lasso, Distributions, ThreadSafeDicts
    # @time include(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")

    # CONSTANTS and trackers
    # --------
    val = :value
    opt["matched"] = 3    # how many matched non-iso cycles to add per iso cycle
    opt["has_df"] = false # tracks state about whether df var for glm is backed up in a jld2 file

    # DATAFRAME-IZE firing rate matrix
    # and add properties of interest
    # --------------------------------
    # Get dataframe of R
    Rdf = DataFrame(R; name=val)
    # Set cycle via start cycle times
    cycles.time = cycles.start
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
    iso_cycles = unique(@subset(Rdf, :isolated_sum .> 0).cycle)
    indexers = [:time,:isolated_sum]
    selector = :area in propertynames(Rdf_cycles) ? Not([:time, :area]) : Not(:time)

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
    @showprogress "annotate cycle with behavior" for cycle in eachrow(cycles)
        inds = DIutils.in_range(beh.time, (cycle.start, cycle.stop))
        for prop in props
            try
            cycle[prop] = nanmean(beh[inds, prop])
            catch
            end
        end
    end
    # Get the grid binning and occupancy of these cycles
    # -----------------------------------------------------
    props = [:x, :y, :speed, :startWell, :stopWell]
    speed_1σ  = Float32(std(beh.speed))
    theta_cycle_req = 1f0 # disable required number of cy
    theta_cycle_dt = median(diff(cycles.start))
    grid_kws =
            (;widths       = [10f0,  10f0,  speed_1σ, 1f0,   1f0],
              radiusinc    = [0.2f0, 0.2f0, 0f0,      0f0,   0f0],
              maxrad       = [6f0,   6f0,   0.4f0,    0.4f0, 0.4f0],
              radiidefault = [2f0,   2f0,   0.4f0,    0.4f0, 0.4f0],
              steplimit=1,
              thresh=theta_cycle_req
             )
    grd = DIutils.binning.get_grid(cycles, props; grid_kws...)
    occ = DIutils.binning.get_occupancy_indexed(Rdf, grd)
    # Save!!!
    # -----------------------------------------------------
    has_df = false
    jldopen(path_iso(opt; append="_cyclewise"), "a") do storage
        storage["grd"] = grd
        storage["occ"] = occ
        has_df = "df" in keys(storage)
    end

    for cyc in iso_cycles

    end
end

function grab_cycle_data(Rdf_cycles::GroupedDataFrame, cyc::Int)::DataFrame
     # Address cycles of interest
     🔑s = [(;cycle=cyc) 
            for cyc in UnitRange(cyc-8, cyc+8)
           ]

    # Grab each cycle of activity
    U = [begin
         u = unstack(Rdf_cycles[🔑], indexers, :unit, val, combine=last) # TODO investigate nonunque
         u = combine(u, selector .=> [mean], renamecols=false)
     end
        for 🔑 in 🔑s if 🔑 in keys(Rdf_cycles)]
     # @info combine(groupby(Rdf_cycles[🔑],:unit),:time=>x->length(x)==length(unique(x)))

    cycs = [🔑.cycle for 🔑 in 🔑s 
                if 🔑 in keys(Rdf_cycles)]
    relcycs = [🔑.cycle-cyc for 🔑 in 🔑s 
                  if 🔑 in keys(Rdf_cycles)]

    # Added df to list
    hcat(DataFrame([cycs,relcycs],[:cycs,:relcycs]), 
                 vcat(U...; cols=:union))
end

# GET CYCLEWISE INFORMATION
fn = path_iso(opt; append="_cyclewise")
if (!isfile(fn) && has_df) || opt["overwrite"]

    df = Vector{Union{Missing,DataFrame}}(missing, length(iso_cycles))
    matched_cycle_holder = Vector{Union{Int,Missing}}(missing, opt["matched"])
    cyc_error =  Dict() 
    Infiltrator.clear_disabled!()
    prog = Progress(length(iso_cycles); desc="cycle df")
    #= Threads.@threads =# for (i,cyc) in collect(enumerate(iso_cycles))
        try
            
            # Push the isolated cycle and its preceding following cycles
            push!(df, grab_cycle_data(Rdf_cycles, cyc))
            matched_cycs = rand(@subset(cycles, :cycle .== cyc).matched,
                                opt["matched"])
            # Push MATCHED cycles
            for mc in matched_cycs
                push!(df, grab_cycle_data(Rdf_cycles, mc))
            end
            next!(prog)
        catch exception
            cyc_error[cyc] = exception
        #     if mod(i, 100) == 0
        #         @info cyc_error
        #     end
            sleep(0.1)
        end
    end
    @info cyc_error
    df = vcat(df[(!).(ismissing.(df))]...)
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

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF SPIKE COUNTS  🔺
# ========================

# Now let's take the formula and apply them
formulae, models   = OrderedDict(), ThreadSafeDict()
formulae["ca1pfc"] = construct_predict_spikecount(df, "CA1");
formulae["pfcca1"] = construct_predict_spikecount(df, "PFC");
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
formulae["ca1pfc"] = construct_predict_iso(df, "CA1", :has);
formulae["pfcca1"] = construct_predict_iso(df, "PFC", :has);
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
formulae["ca1pfc"] = construct_predict_iso(df, "CA1", :count);
formulae["pfcca1"] = construct_predict_iso(df, "PFC", :count);
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
    catch
        models[(;indep, relcyc)] = nothing
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

plot(
    glmplot(model_spikecount, "pfcca1", func, label="spike count, general"),
    glmplot(model_isocount,    "pfcca1", func, label="iso count"),
    glmplot(model_hasiso,     "pfcca1", func, label="iso occurance");
    ylims=(0,1),
    link=:y,
    layout=grid(3,1)
)
