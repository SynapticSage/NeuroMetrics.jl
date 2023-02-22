
# READY DATATYPES AND MODULES

# IMPORTS AND DATA
# --------
using GLM, Lasso, Distributions, ThreadSafeDicts, DrWatson
include(scriptsdir("isolated","load_isolated.jl"))
@time include("../load_isolated.jl")

#   _  _     ____            _        _                     _ _        
# _| || |_  | __ )  __ _ ___(_) ___  (_)___  ___  ___ _ __ (_) | _____ 
#|_  ..  _| |  _ \ / _` / __| |/ __| | / __|/ _ \/ __| '_ \| | |/ / _ \
#|_      _| | |_) | (_| \__ \ | (__  | \__ \ (_) \__ \ |_) | |   <  __/
#  |_||_|   |____/ \__,_|___/_|\___| |_|___/\___/|___/ .__/|_|_|\_\___|
#                                                    |_|               
#     _        _       
# ___| |_ __ _| |_ ___ 
#/ __| __/ _` | __/ __|
#\__ \ || (_| | |_\__ \
#|___/\__\__,_|\__|___/
                      

# Figure out cycles with some number of isolated spikes
begin
    ca1cycstat = combine(groupby(@subset(spikes, :area .== "CA1"), 
                    :cycle), 
    :isolated => sum, 
    :unit => (n->unique(n)) => :diversity,
    [:isolated, :unit] => 
        ((i,u)->(length(unique(u[i.==true])))) => :isodiv,
    [:isolated, :unit] => 
        ((i,u)->(length(unique(u[i.==false])))) => :adjdiv,
    )
    dropmissing!(ca1cycstat, :cycle)
end
begin
    pfccycstat = combine(groupby(@subset(spikes, :area .== "PFC"), 
                    :cycle), 
    :isolated => sum => :pfcisosum, 
    :unit => (n->unique(n)) => :diversity,
    [:isolated, :unit] => 
        ((i,u)->(length(unique(u[i.==true])))) => :pfcisodiv,
    [:isolated, :unit] => 
        ((i,u)->(length(unique(u[i.==false])))) => :pfcadjdiv,
    )
    dropmissing!(pfccycstat, :cycle)
end

cyccellstat = combine(groupby(spikes, [:cycle, :unit]),
    :isolated => sum => :i
)


# Register datasets
@showprogress "registration" for (name,obj) in zip(("cycle","spikes","Rdf"),
                                                    (cycles, spikes, Rdf))
    @info "register" name
    DIutils.filtreg.register(ca1cycstat, obj, on="cycle", 
                     transfer=String.(setdiff(propertynames(ca1cycstat),
                                            [:cycle])))
    DIutils.filtreg.register(pfccycstat, obj, on="cycle", 
                     transfer=String.(setdiff(propertynames(ca1cycstat),
                                            [:cycle])))
end
begin
    css = groupby(cyccellstat, :unit)
    for obj in (spikes, Rdf)
        obj[!,:i] = Vector{Union{Int,Missing}}(missing, size(obj,1))
        tmp = groupby(obj, :unit)
        K = intersect(keys(tmp.keymap), keys(tmp.keymap))
        @showprogress for k in K
            DIutils.filtreg.register(css[k], tmp[k], on="cycle", 
                transfer=["i"])
        end
    end
end
DIutils.filtreg.register(cells, Rdf, on="unit", transfer=["area"])

# Annotate spikes and Rdf
dropmissing!(Rdf, :isolated_sum)
Rdf_cycles    = groupby(Rdf, [:cycle])
indexers = [:time,:isolated_sum]

begin
    kws=(;label="")
    plot(
     (@df ca1cycstat histogram(:isolated_sum;xlabel="iso spikes emitted",kws...)),
     (@df ca1cycstat histogram(:isodiversity;xlabel="# of uniq iso cells",kws...)),
     (@df ca1cycstat histogram(:adjdiversity;ylabel="# of uniq adj cells",kws...)),
     Plot.blank((plot();Plots.annotate!(0.5,0.5,text("Cycle statistics",14))),
            visible=false, size=(100,50)),
     layout=grid(2,2))
end


#    _  _     __  __       _       _       _   _          _        
#  _| || |_  |  \/  | __ _| |_ ___| |__   | |_| |__   ___| |_ __ _ 
# |_  ..  _| | |\/| |/ _` | __/ __| '_ \  | __| '_ \ / _ \ __/ _` |
# |_      _| | |  | | (_| | || (__| | | | | |_| | | |  __/ || (_| |
#   |_||_|   |_|  |_|\__,_|\__\___|_| |_|  \__|_| |_|\___|\__\__,_|
#                                                                  
#                  _           
#   ___ _   _  ___| | ___  ___ 
#  / __| | | |/ __| |/ _ \/ __|
# | (__| |_| | (__| |  __/\__ \
#  \___|\__, |\___|_|\___||___/
#       |___/                  
# 
# ( For matching isolated theta cycles with behavior "equivalent"
#  non-isolated theta cycles)
begin

    # Obtain average of properties of interest per θ cycle
    # -----------------------------------------------------
    for prop in matchprops
        cycles[!,prop] = Vector{Union{Missing,eltype(beh[!,prop])}}(missing,
            size(cycles,1))
    end
    E = Threads.Atomic{Int}(0)
    Threads.@threads for i in collect(1:size(cycles,1))
        cycle = cycles[i,:]
        inds = DIutils.in_range(beh.time, [cycle.start, cycle.stop])
        for prop in matchprops
            try
                cycle[prop] = nanmean(collect(skipmissing(beh[inds, prop])))
            catch
                Threads.atomic_add!(E, 1)
            end
        end
    end
    println("Fail percent ", E[] / size(cycles,1))
    for prop in matchprops
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
    grd = DIutils.binning.get_grid(use_behavior ? beh : cycles, matchprops;
                                    grid_kws...)
    occ = DIutils.binning.get_occupancy_indexed(cycles, grd)
    println("Percent cycles counted ",
            sum(occ.count)/size(dropmissing(cycles,matchprops),1))

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
    iso_cycles = unique(
        @subset(Rdf, :isolated_sum .> 0, :hasocc .== true).cycle)
    cycles.matched = Vector{Union{Vector{Int32}, Missing}}(missing,
        size(cycles,1))
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

#    _  _      ____      _                     _      
#  _| || |_   / ___| ___| |_    ___ _   _  ___| | ___ 
# |_  ..  _| | |  _ / _ \ __|  / __| | | |/ __| |/ _ \
# |_      _| | |_| |  __/ |_  | (__| |_| | (__| |  __/
#   |_||_|    \____|\___|\__|  \___|\__, |\___|_|\___|
#                                   |___/             
#  _           _       _                
# | |__   __ _| |_ ___| |__   ___  ___  
# | '_ \ / _` | __/ __| '_ \ / _ \/ __| 
# | |_) | (_| | || (__| | | |  __/\__ \ 
# |_.__/ \__,_|\__\___|_| |_|\___||___/ 
#
# (each iso/noniso cycle plus precedents and antecedents)
                                      
fn = path_iso(opt; append="_cyclewise")
if (!isfile(fn) && has_df) || opt["overwrite"]

    df = Vector{Any}(undef, 16)
    for thread in 1:Threads.nthreads()
        df[thread] = Vector{Union{Missing,DataFrame}}(missing,
            length(iso_cycles))
    end
    matched_cycle_holder = Vector{Union{Int,Missing}}(missing, opt["matched"])
    cyc_error =  Dict() 
    Infiltrator.clear_disabled!()
    prog = Progress(length(iso_cycles); desc="grabing cycle batches into df")
    Threads.@threads for (i,cyc) in collect(enumerate(iso_cycles))
        try
            tid = Threads.threadid()
            cyc_batch = i
            # Push the isolated cycle and its preceding following cycles
            push!(df[tid], grab_cycle_data(Rdf_cycles, cyc; val, indexers,
                                            cyc_batch, cyc_match=0))
            matched_cycs = @subset(cycles, :cycle .== cyc).matched[1]
            # Push MATCHED cycles
            for (j,mc) in enumerate(matched_cycs)
                push!(df[tid], grab_cycle_data(Rdf_cycles, mc; val, indexers, 
                                                cyc_batch, cyc_match=j))
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
    # TODO not INT because it's gaussian smoothed
    df[:,neuroncols] .*= median(diff(beh.time)) 
    df[:,neuroncols] .= round.(df[:,neuroncols])
    df = transform(df, neuroncols .=> n -> convert(Vector{Int64}, n), 
        renamecols=false)

    # Checkpoint
    rm(fn)
    jldsave(fn; df)

else
    storage = JLD2.jldopen(fn)
    df = storage["df"]
    close(storage)
end

# Clean data frame
col_all_zero = map(v->all(v.==0), eachcol(df))
df = df[!, Not(names(df)[col_all_zero])]


function get_dx_dy(df, relcyc)
    dx = @subset(df, :relcycs .== relcyc)
    dy = @subset(df, :relcycs .== 0)
    dx = groupby(dx, [:cyc_batch, :cyc_match])
    dy = groupby(dy, [:cyc_batch, :cyc_match])
    kx, ky = keys(dx.keymap), keys(dy.keymap)
    k = intersect(kx,ky)
    dx = sort(vcat([dx[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    dy = sort(vcat([dy[kk] for kk in k]...), [:cyc_batch, :cyc_match])
    dx, dy
end
@time dx_dy = Dict(relcyc => get_dx_dy(df, relcyc) for relcyc in -8:8)

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF SPIKE COUNTS  🔺
# ========================

# # Now let's take the formula and apply them
# formulae, models, cache   = OrderedDict(), ThreadSafeDict(), ThreadSafeDict()
# formulae["ca1pfc"] = construct_predict_spikecount(df, cells, "CA1");
# formulae["pfcca1"] = construct_predict_spikecount(df, cells, "PFC");
# @assert !isempty(first(values(formulae)))
# glmsets = []
# for indep in ("ca1pfc", "pfcca1"), relcyc in -8:8, f in formulae[indep]
#     push!(glmsets, (indep, relcyc, f))
# end
#
# using GLMNet, MultivariateStats, MLJ, ScikitLearn
# using MATLAB
#
# prog = Progress(length(glmsets); desc="GLM spike counts")
# Threads.@threads for (indep, relcyc, f) in glmsets
#     unit = parse(Int,string(f.lhs))
#     try
#         dx, dy = dx_dy[relcyc]
#         _, XX = modelcols(apply_schema(f, schema(f, dx)), dx)
#         y, _  = modelcols(apply_schema(f, schema(f, dy)), dy)
#         FittedPoisson = fit_mle(Poisson, Int.(y))
#         # TODO better research glmnet -- do i need to link manually -- its
#         # negative output for a neural firing prediction
#         models[(;indep, relcyc, unit)] = m =  glm(XX, y, FittedPoisson)
#         # models[(;indep, relcyc, unit)] = m =  glmnetcv(XX, y,FittedPoisson)
#         #  mat"[$B, $stats] = lassoglm(double($XX), double($y), 'poisson', 
#         #   'alpha', 0.5, 'CV', 3, 'MCReps', 5, 'link', 'log')"
#         # models[(;indep, relcyc)] = (;B, stats)
#         cache[(;indep, relcyc, unit)] = (;XX, y)
#     catch exception
#         models[(;indep, relcyc, unit)] = exception
#         sleep(0.1)
#     end
#     next!(prog)
# end
#
# @info "saving spikcount"
# if isdefined(Main, :storage); close(storage); end
# storage = jldopen(fn, "a")
# "model_spikecount" in keys(storage) ? delete!(storage, "model_spikecount")
#           : nothing
# global storage["model_spikecount"] = model_spikecount = models
# @info "saved"
# if isdefined(Main, :storage); close(storage); end

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
        dx, dy = dx_dy[relcyc]
        _, XX = modelcols(apply_schema(f, schema(f, dx)), dx)
        y, _  = modelcols(apply_schema(f, schema(f, dy)), dy)
        FittedBinomial = fit_mle(Binomial, 1, y)
        models[(;indep, relcyc)] = m =  glm(XX, y, FittedBinomial)
        # models[(;indep, relcyc)] = glm(XX, y, FittedBinomial)
        # mat"[$B, $stats] = lassoglm(double($XX), double($y), 'binomial',
        # 'alpha', 0.5, 'CV', 3, 'MCReps', 5, 'link', 'logit')"
        # models[(;indep, relcyc)] = (;B, stats)
    catch exception
        models[(;indep, relcyc)] = exception
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
        dx, dy = dx_dy[relcyc]
        _, XX = modelcols(apply_schema(f, schema(f, dx)), dx)
        y, _  = modelcols(apply_schema(f, schema(f, dy)), dy)
        FittedPoisson = fit_mle(Poisson, Int.(y))
        models[(;indep, relcyc)] = m =  fit(GLM.GeneralizedLinearModel, XX, y,
                                            FittedPoisson)
        # mat"[$B, $stats] = lassoglm(double($XX), double($y), 'binomial',
        # 'alpha', 0.5, 'CV', 3, 'MCReps', 5, 'link', 'logit')"
        # models[(;indep, relcyc)] = (;B, stats)
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


#    _  _     ____  _       _   _   _             
#  _| || |_  |  _ \| | ___ | |_| |_(_)_ __   __ _ 
# |_  ..  _| | |_) | |/ _ \| __| __| | '_ \ / _` |
# |_      _| |  __/| | (_) | |_| |_| | | | | (_| |
#   |_||_|   |_|   |_|\___/ \__|\__|_|_| |_|\__, |
#                                           |___/ 
#
Plot.setfolder("isolated", "glm")
Plot.setappend((;animal, day, tet))

using Metrics, .Table, GLMNet
import Dates
func = (k,v)->if v isa Exception || v === nothing
    NaN
else
    # adjr2(v, :devianceratio)
    if v isa GLMNet.GLMNetCrossValidation
        XX, y = cache[k]
        y0 = GLMNet.predict(v, XX)
        #Metrics.r2_score(y0, y)
        Dict("y"=>y, "y0"=>y0)
    elseif v isa StatisticalModel
        adjr2(v, :devianceratio)
    else
        NaN
    end
end
D = to_dataframe(model_spikecount, func)

Dsum = sort(combine(groupby(D, [:relcyc, :indep]), :value=>nanmedian=>:value),
        [:indep,:relcyc])
@df @subset(Dsum,:indep .== "pfcca1") begin
    scatter(:relcyc, :value, alpha=0.5, ylim=(0,1))
end

P=[(@df d scatter(:relcyc, :value, alpha=0.5, ylim=(0,1)))
    for d in groupby(D, :unit)]

plot(
    plot(
        glmplot(model_isocount, "pfcca1", func, label="iso count"),
        glmplot(model_hasiso,   "pfcca1", func, label="iso occurance");
        link=:y,
        layout=grid(2,1),
        title= "pfc -> ca1"
    ),
    plot(
        glmplot(model_isocount, "ca1pfc", func, label="iso count"),
        glmplot(model_hasiso,   "ca1pfc", func, label="iso occurance");
        link=:y,
        layout=grid(2,1),
        title="ca1 -> pfc"
    ),
    ylims=(0,0.4),
    size=(600, 900)
)
Plot.save("iso, count and occur, date=$(Dates.now())")
