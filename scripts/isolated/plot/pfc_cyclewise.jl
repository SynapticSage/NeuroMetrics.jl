# IMPORTS AND DATA
# --------
checkpoint = false
if checkpoint
    include(scriptsdir("isolated","load_cyclewise_checkpoint.jl"))
    @time include("../load_cyclewise_checkpoint.jl")
else
    @time include(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")
end
include("./imports_isolated.jl")
include("../imports_isolated.jl")


# Filtration
pyramidal = cells[cells.meanrate .< 5, :unit]
spikecount = combine(groupby(spikes, :unit), nrow=>:spikecount, 
                             :isolated=>sum=>:isospikecount)
spikecount_criteria = 
    spikecount[spikecount.isospikecount .> 5 .&& 
        spikecount.spikecount .> 50,:unit]
neurons_of_interest = intersect(pyramidal, spikecount_criteria)
spikes_sub = @subset(spikes, :unit .∈ (neurons_of_interest,))
Rdf_sub    = @subset(Rdf, :unit .∈ (neurons_of_interest,))
cells_sub  = @subset(cells, :unit .∈ (neurons_of_interest,))
@assert size(spikes,1) != size(spikes_sub,1)

#   _  _     ____            _        _                     _ _        
# _| || |_  | __ )  __ _ ___(_) ___  (_)___  ___  ___ _ __ (_) | _____ 
#|_  ..  _| |  _ \ / _` / __| |/ __| | / __|/ _ \/ __| '_ \| | |/ / _ \
#|_      _| | |_) | (_| \__ \ | (__  | \__ \ (_) \__ \ |_) | |   <  __/
#  |_||_|   |____/ \__,_|___/_|\___| |_|___/\___/|___/ .__/|_|_|\_\___|
# ___| |_ __ _| |_ ___ 
#/ __| __/ _` | __/ __|
#\__ \ || (_| | |_\__ \
#|___/\__\__,_|\__|___/
                      

# Figure out cycles with some number of isolated spikes
begin
    ca1cycstat = combine(groupby(@subset(spikes_sub, :area .== "CA1"), 
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
    pfccycstat = combine(groupby(@subset(spikes_sub, :area .== "PFC"), 
                    :cycle), 
    :isolated => sum => :pfcisosum, 
    :unit => (n->unique(n)) => :pfcdiversity,
    [:isolated, :unit] => 
        ((i,u)->(length(unique(u[i.==true])))) => :pfcisodiv,
    [:isolated, :unit] => 
        ((i,u)->(length(unique(u[i.==false])))) => :pfcadjdiv,
    )
    dropmissing!(pfccycstat, :cycle)
end

cyccellstat = combine(groupby(spikes_sub, [:cycle, :unit]),
    :isolated => sum => :i
)


# Register datasets
@showprogress "registration" for (name,obj) in zip(("cycle","spikes","Rdf"),
                                                    (cycles, spikes_sub, Rdf_sub))
    @info "register" name
    DIutils.filtreg.register(ca1cycstat, obj, on="cycle", 
                     transfer=String.(setdiff(propertynames(ca1cycstat),
                                            [:cycle])))
    DIutils.filtreg.register(pfccycstat, obj, on="cycle", 
                     transfer=String.(setdiff(propertynames(pfccycstat),
                                            [:cycle])))
end
begin
    css = groupby(cyccellstat, :unit)
    for obj in (spikes_sub, Rdf_sub)
        obj[!,:i] = zeros(Union{Missing,Int64}, size(obj,1))
        tmp = groupby(obj, :unit)
        K = intersect(keys(tmp.keymap), keys(tmp.keymap))
        @showprogress for k in K
            inds_css = (!).(ismissing.(css[k][!, :cycle]))
            inds_tmp = (!).(ismissing.(tmp[k][!, :cycle]))
            a,b = DIutils.filtreg.register(css[k][inds_css,:],
                convert_type=Int64,
                tmp[k][inds_tmp,:],  on="cycle", transfer=["i"])
            css[k][inds_css,:], tmp[k][inds_tmp,:] = a, b
        end
    end
end
DIutils.filtreg.register(cells, Rdf_sub, on="unit", transfer=["area"])

# Annotate spikes and Rdf
dropmissing!(Rdf_sub, :isolated_sum)
iso_cycles = unique(@subset(Rdf_sub, :isolated_sum .> 0, 
                    :hasocc .== true).cycle)
indexers = [:time, :isolated_sum, :pfcisosum, :bins]
commit_vars()

# 
#    _  _     ____                                
#  _| || |_  |  _ \ _ __ ___ _ __   __ _ _ __ ___ 
# |_  ..  _| | |_) | '__/ _ \ '_ \ / _` | '__/ _ \
# |_      _| |  __/| | |  __/ |_) | (_| | | |  __/
#   |_||_|   |_|   |_|  \___| .__/ \__,_|_|  \___|
#                           |_|                   
#  _ __ ___   __ _| |_ ___| |__   | |_| |__   ___| |_ __ _ 
# | '_ ` _ \ / _` | __/ __| '_ \  | __| '_ \ / _ \ __/ _` |
# | | | | | | (_| | || (__| | | | | |_| | | |  __/ || (_| |
# |_| |_| |_|\__,_|\__\___|_| |_|  \__|_| |_|\___|\__\__,_|
#                                                          
# ( For matching isolated theta cycles with behavior "equivalent"
#  non-isolated theta cycles)
begin

    # (1) Place behavior into cycles
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

    # (2) Get an adaptive grid
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
end

# Match cycles
match_cycles!(cycles, Rdf_sub, occ; matches=opt["matched"], iso_cycles)

# Save!!! cyclewise specific
# -----------------------------------------------------
commit_cycwise_vars()

#    _  _      ____      _                     _      
#  _| || |_   / ___| ___| |_    ___ _   _  ___| | ___ 
# |_  ..  _| | |  _ / _ \ __|  / __| | | |/ __| |/ _ \
# |_      _| | |_| |  __/ |_  | (__| |_| | (__| |  __/
#   |_||_|    \____|\___|\__|  \___|\__, |\___|_|\___|
# | |__   __ _| |_ ___| |__   ___  ___ | |
# | '_ \ / _` | __/ __| '_ \ / _ \/ __| _|
# | |_) | (_| | || (__| | | |  __/\__ \ 
# |_.__/ \__,_|\__\___|_| |_|\___||___/ 
#

# Prepare for a 10 minute bucketed approach
dT = diff(beh.time)
bins = Int(maximum(floor.(cumsum(dT)/60/10))) #get number of 10 minute buckets
beh.bins = DIutils.binning.digitize(beh.time, bins);
checkbins(:beh)
DIutils.filtreg.register(beh, cycles, on="time",
    transfer=["epoch","bins"]);
checkbins(:cycles);
df.cycle = df.cyc_central;
# DIutils.filtreg.register(cycles, df, on="cycle",
#     transfer=["epoch","time","bins"]);
DIutils.filtreg.register(cycles, Rdf_sub, on="cycle", transfer=["epoch","bins"]);
checkbins(:Rdf_sub);

# (each iso/noniso cycle plus precedents and antecedents)
Rdf_cycles  = groupby(Rdf_sub, [:cycle])
df, cyc_errors = df_FRpercycle_and_matched(cycles, Rdf_cycles,
    beh, val; iso_cycles=iso_cycles, indexers, cycrange=opt["cycles"])
checkbins(:df);
# Checkpoint
commit_cycwise_vars()

#    _  _      ____           _          
# _| || |_   / ___|__ _  ___| |__   ___ 
#|_  ..  _| | |   / _` |/ __| '_ \ / _ \
#|_      _| | |__| (_| | (__| | | |  __/
#  |_||_|    \____\__,_|\___|_| |_|\___|
#                                       
dx_dy = ThreadSafeDict()
Threads.@threads for relcyc in collect(-opt["cycles"]:opt["cycles"])
    dx_dy[(;relcyc)] = get_dx_dy(df, relcyc)
end
dx_dy = Dict(dx_dy)

DFB= groupby(df, :bins)
prog = Progress((opt["cycles"]*2+1) * length(DFB))
@time dx_dy_bin = ThreadSafeDict()
sets = [(relcyc, DFB[b], b) 
    for relcyc in -opt["cycles"]:opt["cycles"],
    b in eachindex(DFB)]
Threads.@threads for (relcyc, dfb, b) in sets
    key = (;relcyc, bin=b[1])
    try
        dx_dy_bin[key] = get_dx_dy(dfb, relcyc)
    catch exception
        @warn "$key failed" exception
    end
    next!(prog)
end
dx_dy_bin = Dict(dx_dy_bin)

# @time dx_dy = merge(dx_dy, get_futurepast_blocks(df))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
include(scriptsdir("isolated","imports_isolated.jl"))

for col in eachcol(df[!,[x for x in names(df) if occursin("_i", x)]] )
    replace!(col, missing=>0)
    col = convert(Vector{Union{Missing,Float64}}, col)
end
df[!,[x for x in names(df) if occursin("_i", x)]]

# =========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF HAS isolated spike per cell  🔺
# ========================
model_cellhasiso, cacheiso, shuffle_cellhasiso = 
                            initorget("model_cellhasiso"), ThreadSafeDict(),
                            initorget("shuffle_cellhasiso"; obj=Dict())
model_cellhasiso = OrderedDict()
pos = (df, dx_dy, cells, construct_predict_isospikecount, :pyglm)
pos = (df, Dict(first(dx_dy)), cells, construct_predict_isospikecount, :pyglm)
kws = (Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""],
    ytrans=x->Float64(x>0), xtrans=x->Float64(x), modelz=model_cellhasiso,
    )
# isolated.run_glm!(pos...;kws...) 

# Relative cycle 0 runs of the model
isolated.run_glm!(df, dx_dy, cells, construct_predict_spikecount, :pyglm;
    Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""],
    ytrans=x->Float64(x>0), xtrans=x->Float64(x), modelz=model_cellhasiso,
    )

#    _  _     _____                _     _           _                 _   
#  _| || |_  |_   _| __ ___  _   _| |__ | | ___  ___| |__   ___   ___ | |_ 
# |_  ..  _|   | || '__/ _ \| | | | '_ \| |/ _ \/ __| '_ \ / _ \ / _ \| __|
# |_      _|   | || | | (_) | |_| | |_) | |  __/\__ \ | | | (_) | (_) | |_ 
#   |_||_|     |_||_|  \___/ \__,_|_.__/|_|\___||___/_| |_|\___/ \___/ \__|
# TESTING WITH THE ENTIRE DATAFRAME
model = Dict()

# Dry run of the model
isolated.run_glm!(df, dx_dy, cells, construct_predict_spikecount, :matlab;
    Dist=Distributions.Poisson(), unitwise=true, 
    xtrans=x->Float64(x), ytrans=x->Float64(x),
    modelz=model, dryrun=true)

# Test on indvidual formulae
fca1pfc = construct_predict_spikecount(df, cells, "PFC")
fpfcca1 = construct_predict_spikecount(df, cells, "CA1")
glm_list_df = []
@showprogress "df, run types" for type in [:matlab, :pyglm, :r]
    A = glm_(fca1pfc, df, df, Poisson(); type, handle_exception=:warn,
        desc=Dict("desc"=>"CA1->PFC", "typ"=>type))
    push!(glm_list_df, A)
    B = glm_(fpfcca1, df, df, Poisson(); type, handle_exception=:warn,
        desc=Dict("desc"=>"PFC->CA1", "typ"=>type))
    push!(glm_list_df, B)
end

# # Running for entire df (not by relcyc)
# isolated.run_glm!(df, cells, construct_predict_spikecount, :pyglm;
#     Dist=Distributions.Binomial(), unitwise=true, modelz=model)

# Testing full firing rate matrix, sans cycles
Rdf_unstack = unstack(Rdf_sub, :time, :unit, :value, combine=mean)
# Test on indvidual formulae
fca1pfc = construct_predict_spikecount(df, cells, "PFC")
fpfcca1 = construct_predict_spikecount(df, cells, "CA1")
glm_list_rdf = []
@showprogress "Rdf, run types" for type in [:matlab, :pyglm, :r]
    A = glm_(fca1pfc, Rdf_unstack, Rdf_unstack, Poisson(); type, 
        handle_exception=:warn, desc=Dict("desc"=>"CA1->PFC", "typ"=>type))
    push!(glm_list_rdf, A)
    B = glm_(fpfcca1, Rdf_unstack, Rdf_unstack, Poisson(); type, 
        handle_exception=:warn, desc=Dict("desc"=>"PFC->CA1", "typ"=>type))
    push!(glm_list_rdf, B)
end

#    _  _     _____                  _ _ _         
#  _| || |_  | ____|__ _ _   _  __ _| (_) |_ _   _ 
# |_  ..  _| |  _| / _` | | | |/ _` | | | __| | | |
# |_      _| | |__| (_| | |_| | (_| | | | |_| |_| |
#   |_||_|   |_____\__, |\__,_|\__,_|_|_|\__|\__, |
#                     |_|                    |___/ 
#   ___| |__   ___  ___| | _____ 
#  / __| '_ \ / _ \/ __| |/ / __|
# | (__| | | |  __/ (__|   <\__ \
#  \___|_| |_|\___|\___|_|\_\___/

# Comparing the results from glming on Rdf_unst1ack and df
    # Convert res to dataframe
    glm_df = Table.to_dataframe(glm_list_df)
    glm_Rdf = Table.to_dataframe(glm_list_rdf)
    # Plot and compare the R2 values using
#
    # R2 value equality?
    H = HypothesisTests.KruskalWallisTest(glm_Rdf[:,:R2], glm_df[:,:R2], 0.05)
    @df glm_Rdf histogram(:R2, group=:desc, label=:desc, 
                            title="R2 values for Rdf")
    @df glm_df histogram(:R2, group=:desc, label=:desc,
                            title="R2 values for df")

    # Coefficient equality?
    H = HypothesisTests.KruskalWallisTest(glm_Rdf[:,:coef], glm_df[:,:coef], 0.05)
    @df glm_Rdf histogram(:coef, group=:desc, label=:desc, 
                            title="Coef values for Rdf")
    @df glm_df histogram(:coef, group=:desc, label=:desc,
                            title="Coef values for df")


# tmp = shuffle_cellhasiso
# shuffle_cellhasiso = ThreadSafeDict()
# for (k,v) in tmp
#     push!(shuffle_cellhasiso,k=>v)
# end


# -----------------------------------------------
# Description of the different construct_ methods
# -----------------------------------------------
#
# construct_predict_isospikecount :: 
#     Predicts the number of isolated spikes in a cycle
# Type :: Poisson
#
# construct_predict_iso :: 
#     Whether there is an isolated spike in a cycle
# Type :: Binomial
#
# construct_predict_spikecount :: 
#     Predicts the number of spikes in a cycle
# Type :: Poisson
#
# -----------------------------------------------
# Description of different model variables stored
# in the local JLD2 file
# -----------------------------------------------
#
# model_cellhasiso :: 
#     Model for predicting whether there is an isolated spike in a cycle
#
# model_cellhasiso_matlab ::
#     Model for predicting whether there is an isolated spike in a cycle
#     using the matlab implementation
#
# 
# shuffle_cellhasiso :: 
#     Shuffled model for predicting whether there is an isolated spike in a cycle
# 
# model_spikecount :: 
#     Model for predicting the number of spikes in a cycle
#
# model_counthasiso :: 
#     Model for predicting whether there is an isolated spike in a cycle
#     using the number of spikes in a cycle as a predictor
# 
#

run_shuffle!(pos...; kws..., 
    shuffle_models=shuffle_cellhasiso,
    shufcount=30)

commit_cycwise_vars()
