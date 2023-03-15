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
include(scriptsdir("isolated","imports_isolated.jl"))
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



    cycles.hasocc = (!).(ismissing.(occ.datainds))
    DIutils.filtreg.register(cycles, Rdf, on="cycle", transfer=["hasocc"])
    DIutils.filtreg.register(cycles, Rdf_sub, on="cycle", transfer=["hasocc"])
    iso_cycles = unique(@subset(Rdf_sub, :isolated_sum .> 0, 
                        :hasocc .== true).cycle)
    indexers = [:time, :isolated_sum, :pfcisosum, :bins]

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
# DIutils.filtreg.register(cycles, df, on="cycle",
#     transfer=["epoch","time","bins"]);
DIutils.filtreg.register(cycles, Rdf_sub, on="cycle", 
    transfer=["epoch","bins"]);
DIutils.filtreg.register(cycles, spikes, on="cycle", 
    transfer=["bins"]);
checkbins(:Rdf_sub);

# (each iso/noniso cycle plus precedents and antecedents)
Rdf_cycles  = groupby(Rdf_sub, [:cycle])
df, cyc_errors = df_FRpercycle_and_matched(cycles, Rdf_cycles,
    beh, val; iso_cycles=iso_cycles, indexers, cycrange=opt["cycles"])
df.cycle = df.cyc_central;
checkbins(:df);

# Remove missing vals from isolated spike count columns
for col in eachcol(df[!,[x for x in names(df) if occursin("_i", x)]] )
    replace!(col, missing=>0)
    col = convert(Vector{Union{Missing,Float64}}, col)
end

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sets to iterate GLM over
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_spikecount = Dict()
model_cellhasiso = Dict()
B = [b for b in collect(keys(dx_dy_bin))]

uAreas = unique(cells.area)
ind_dep = [(ind, dep) for dep in uAreas, ind in uAreas]
bf(F) = [(b, f) for b in B, f in F]
begin # test first sets
    (dep, ind) = first(ind_dep)
    (b, f) = first(bf(construct_predict_spikecount(df, cells, ind;
                dep_area=dep)))
    (b, f) = first(bf(construct_predict_isospikecount(df, cells, ind;
                dep_area=dep)))
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # =============================================
    # PREDICT SPIKCOUTNS ALL TIMES BY 10 minute BIN
    # =============================================
    @showprogress "area-interact" for (dep, ind) in ind_dep
        F = construct_predict_spikecount(df, cells, ind; dep_area=dep)
    @showprogress "single run" for (b,f) in bf(F)
        XXb, yb = dx_dy_bin[b]
        key = (;b..., unit=parse(Int,string(f.lhs)), dir="$ind => $dep")
        if !opt["overwrite"] && haskey(model_spikecount, key)
            continue
        end
        m = glm_(f, XXb, yb, "poisson"; type=:pyglm,
            desc=Dict("bin"=>key.bin, "unit"=>key.unit, "relcyc"=>key.relcyc,
                      "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
                "dep"=>dep, "ind"=>ind))
        model_spikecount[key] = m
    end
    end
    commit_cycwise_vars("model_spikecount")

    shuffle_spikecount = initorget("shuffle_spikecount")
    @showprogress "shuffle" for _ in 1:opt["shuffle"]
        shuf = hash(rand())
    @showprogres "area-interact" for (dep, ind) in ind_dep
        F = construct_predict_spikecount(df, cells, ind; dep_area=dep)
    @showprogress "ones shuf" for (b,f) in bf(F)
        XXb, yb = dx_dy_bin[b]
        key = (;shuf, relcyc=0, bin=b.bin, unit=parse(Int,string(f.lhs)))
        if !opt["overwrite"] && haskey(shuffle_spikecount, key)
            continue
        end
        s = shuffle_glm_(f, XXb, yb, "poisson"; type=:pyglm,
            desc=Dict("bin"=>key.bin, "unit"=>key.unit, "relcyc"=>key.relcyc,
                      "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
                "dep"=>dep, "ind"=>ind))
        shuffle_spikecount[key] = m
    end
    end
    end
    commit_cycwise_vars("shuffle_spikecount")

    # =============================================
    # PREDICT CELL HAS ISO SPIKE BY 10 minute BIN
    # =============================================
    @showprogress "area-interact" for (dep, ind) in ind_dep
        F = construct_predict_isospikecount(df, cells, ind; dep_area=dep)
    @showprogress "single run" for (b,f) in bf(F)
        XXb, yb = dx_dy_bin[b]
        nonmiss = vec(any(
                 (!).(ismissing.(Array(yb)) .|| ismissing.(Array(XXb)) ),
            dims=2))
        XXb, yb = XXb[nonmiss, :], yb[nonmiss,:]
        yb[:,uicellcols(yb)] = Int.(yb[!,uicellcols(yb)] .> 0)
        unit=parse(Int, replace(string(f.lhs), "_i"=>""))
        b = (;b..., unit, dir="$ind => $dep", dist="binomial")
        if !opt["overwrite"] && haskey(model_cellhasiso, b)
            continue
        end
        m = glm_(f, XXb, yb, "binomial"; type=:pyglm,
            desc=Dict("bin"=>b.bin, "unit"=>b.unit, "relcyc"=>b.relcyc,
                      "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
                "dep"=>dep, "ind"=>ind))
        model_cellhasiso[b] = m
        sleep(0.001)
    end
    end
    commit_cycwise_vars("model_cellhasiso")

df_iso = Table.to_dataframe(model_cellhasiso, name="value", key_name="prop")
commit_cycwise_vars("df_iso")

@df @subset(df_iso, :prop .== "pseudo_R2") begin
    histogram(:value, bins=100, alpha=0.5, label="pseudo_R2",
        title="positive fraction = $(mean(:value .> 0))")
end
# BINS
@df @subset(df_iso, :prop .== "pseudo_R2") begin
    scatter(:relcyc, :value, alpha=0.1, label="pseudo_R2",
        title="positive fraction = $(mean(:value .> 0))")
end
# BINS
@df @subset(df_iso, :prop .== "pseudo_R2") begin
    scatter(:bin, :value, alpha=0.1, label="pseudo_R2",
        title="positive fraction = $(mean(:value .> 0))")
end
# BINS WITH POSITIVE PSEUDO R2
d = @subset(df_iso, :prop .== "pseudo_R2", :value .> 0.0)
scatter(d.bin, d.value, group=d.unit, alpha=0.5, label="",
    title="positive fraction = $(mean(:value .> 0))",
    ylims=(-0.001, 0.1))
t = Plots.text.(d.dir, [8], ["black"])
p=annotate!.(d.bin, d.value, t)


# =========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF HAS isolated spike per cell  🔺
# ========================
# model_cellhasiso, cacheiso, shuffle_cellhasiso = 
#                             initorget("model_cellhasiso"), ThreadSafeDict(),
#                             initorget("shuffle_cellhasiso"; obj=Dict())
# model_cellhasiso = OrderedDict()
# pos = (df, dx_dy, cells, construct_predict_isospikecount, :pyglm)
# pos = (df, Dict(first(dx_dy)), cells, construct_predict_isospikecount, :pyglm)
# kws = (Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""],
#     ytrans=x->Float64(x>0), xtrans=x->Float64(x), modelz=model_cellhasiso,
#     )
# isolated.run_glm!(pos...;kws...) 

# Relative cycle 0 runs of the model
# isolated.run_glm!(df, dx_dy, cells, construct_predict_spikecount, :pyglm;
#     Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""],
#     ytrans=x->Float64(x>0), xtrans=x->Float64(x), modelz=model_cellhasiso,
#     )


# # Running for entire df (not by relcyc)
# isolated.run_glm!(df, cells, construct_predict_spikecount, :pyglm;
#     Dist=Distributions.Binomial(), unitwise=true, modelz=model)


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
