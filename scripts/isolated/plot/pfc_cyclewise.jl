using GoalFetchAnalysis.Field.metrics: skipnan
# IMPORTS AND DATA
# --------
checkpoint = true
if checkpoint
    include(scriptsdir("isolated","load_cyclewise_checkpoint.jl"))
    @time include("../load_cyclewise_checkpoint.jl")
else
    include(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")
end
include("../imports_isolated.jl")

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
    :unit => (n->unique(n)) => :pfcdiversity,
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
                     transfer=String.(setdiff(propertynames(pfccycstat),
                                            [:cycle])))
end
begin
    css = groupby(cyccellstat, :unit)
    for obj in (spikes, Rdf)
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
DIutils.filtreg.register(cells, Rdf, on="unit", transfer=["area"])

# Annotate spikes and Rdf
dropmissing!(Rdf, :isolated_sum)
Rdf_cycles    = groupby(Rdf, [:cycle])
indexers = [:time, :isolated_sum, :pfcisosum]
commit_vars()

begin
kws=(;label="")

p1=plot(
 (@df ca1cycstat histogram(:isolated_sum;xlabel="iso spikes emitted",kws...)),
 (@df ca1cycstat histogram(:isodiv;xlabel="# of uniq iso cells",kws...)),
 (@df ca1cycstat histogram(:adjdiv;ylabel="# of uniq adj cells",kws...)),
 Plot.blank((plot();Plots.annotate!(0.5,0.5,text("CA1 cycle\nstatistics",14))),
        visible=false, size=(100,50)),
 layout=grid(2,2));

p2=plot(
 (@df pfccycstat histogram(:pfcisosum;xlabel="iso spikes emitted",kws...)),
 (@df pfccycstat histogram(:pfcisodiv;xlabel="# of uniq iso cells",kws...)),
 (@df pfccycstat histogram(:pfcadjdiv;ylabel="# of uniq adj cells",kws...)),
 Plot.blank((plot();Plots.annotate!(0.5,0.5,text("PFC cycle\nstatistics",14))),
            visible=false, size=(100,50)),
    layout=grid(2,2)
);

    plot(p1,p2, size=(1000,500))

end


#    _  _     __  __       _       _       _   _          _        
#  _| || |_  |  \/  | __ _| |_ ___| |__   | |_| |__   ___| |_ __ _ 
# |_  ..  _| | |\/| |/ _` | __/ __| '_ \  | __| '_ \ / _ \ __/ _` |
# |_      _| | |  | | (_| | || (__| | | | | |_| | | |  __/ || (_| |
#   |_||_|   |_|  |_|\__,_|\__\___|_| |_|  \__|_| |_|\___|\__\__,_|
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

    # Save!!! cyclewise specific
    # -----------------------------------------------------
    commit_cycwise_vars()

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

    # Save!!! isolated general vars again after we've added info
    # -----------------------------------------------------------
    commit_vars()

end

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
# (each iso/noniso cycle plus precedents and antecedents)
fn = path_iso(opt; append="_cyclewise")
if (!isfile(fn) && has_df) || opt["overwrite"]

    threading = true
    if threading
        nthreads = Threads.nthreads()
    else
        nthreads = 1
    end
    df = Vector{Vector{Union{Missing,DataFrame}}}(undef, nthreads)
    for thread in 1:nthreads
        df[thread] = Vector{Union{Missing,DataFrame}}()
    end
    matched_cycle_holder = Vector{Union{Int,Missing}}(missing, opt["matched"])
    cyc_error =  Dict() 
    Infiltrator.clear_disabled!()
    prog = Progress(length(iso_cycles); desc="grabing cycle batches into df")
    V = [val, :i]
    E, M = Threads.Atomic{Int}(0), Threads.Atomic{Int}(0)
    #[(length(iso_cycles)-100):end]
    Threads.@threads for (i,cyc) in collect(enumerate(iso_cycles))
    unit = parse(Int,replace(string(f.lhs), "_i"=>""))
        try
            tid = threading ? Threads.threadid() : 1
            cyc_batch = i
            # Push the isolated cycle and its preceding following cycles
            push!(df[tid], grab_cycle_data(Rdf_cycles, cyc, V; indexers,
                                            cycrange=opt["cycles"],
                                            cyc_batch, cyc_match=0))
            matched_cycs = @subset(cycles, :cycle .== cyc).matched[1]
            # Push MATCHED cycles
            if isempty(matched_cycs)
                Threads.atomic_add!(M, 1)
                continue
            end
            for (j,mc) in enumerate(matched_cycs)
                push!(df[tid], 
                    grab_cycle_data(Rdf_cycles, mc, V; indexers, 
                                    cycrange=opt["cycles"],
                                    cyc_batch, cyc_match=j))
            end
            next!(prog)
        catch exception
            cyc_error[cyc] = exception
            Threads.atomic_add!(E, 1)
        #     if mod(i, 100) == 0
        #         @info cyc_error
        #     end
            sleep(0.05)
            next!(prog)
        end
    end
    printstyled("Cycles without match ", M[]/length(iso_cycles), 
          "\nErrored cycles ", E[]/length(iso_cycles), color=:blink)
    vcatnonmiss(df) = vcat(df[(!).(ismissing.(df))]...)
    df = vcatnonmiss.(df)
    df = vcatnonmiss(df)
    @assert :cyc_match ∈ propertynames(df) ||
        unique(df.cyc_match)>1 "FUCK"
    df.has_iso = df.isolated_sum .> 0
    # Spike count
    neuroncols = names(df)[tryparse.(Int, names(df)) .!== nothing]
    # TODO not INT because it's gaussian smoothed
    df[:,neuroncols] .*= median(diff(beh.time)) 
    df[:,neuroncols] .= round.(df[:,neuroncols])
    df = transform(df, neuroncols .=> n -> convert(Vector{Int64}, n), 
        renamecols=false)
    # Clean data frame
    col_all_zero = map(v->all(skipmissing(v.==0)), eachcol(df))
    df = df[!, Not(names(df)[col_all_zero])]
    # Checkpoint
    commit_cycwise_vars()

else
    storage = JLD2.jldopen(fn)
    df = storage["df"]
    close(storage)
end

function get_dx_dy(df, relcyc::Int)
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
function get_futurepast_blocks(df)
    dxf = @subset(df, :relcycs .> 0)
    dxp = @subset(df, :relcycs .<= 0)
    dy = @subset(df, :relcycs .== 0)
    dxp = groupby(dxp, [:cyc_batch, :cyc_match])
    dxf = groupby(dxf, [:cyc_batch, :cyc_match])
    dy = groupby(dy, [:cyc_batch, :cyc_match])
    kxp, kxf, ky = keys(dxp.keymap), keys(dxf.keymap), keys(dy.keymap)
    k = intersect(kxf,kxp,ky)
    dxf = sort(vcat([dxf[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    dxp = sort(vcat([dxp[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    dy  = sort(vcat([dy[kk] for kk in k]...), [:cyc_batch, :cyc_match], )
    Dict("pastcurr"=>(dxp, dy), "future"=>(dxf, dy))
end
@time dx_dy = Dict(relcyc => get_dx_dy(df, relcyc) for relcyc in 
    -opt["cycles"]:opt["cycles"])
@time dx_dy = merge(dx_dy,
    get_futurepast_blocks(df))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
include(scriptsdir("isolated","imports_isolated.jl"))


for col in eachcol(df[!,[x for x in names(df) if occursin("_i", x)]] )
    replace!(col, missing=>0)
    col = convert(Vector{Union{Missing,Float64}}, col)
end
df[!,[x for x in names(df) if occursin("_i", x)]]

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF SPIKE COUNTS  🔺
# ========================
"""
Shortcut function that handles applying my statsmodels formulae to each
cut of the data andn running GLM with the chosen `glmtool`
"""
function run_glm!(runname::String, 
    formula_method::Function, glmtool; Dist=Binomial(),
    unitwise::Bool, unitrep=[], xtrans::Function=identity,
    ytrans::Function=identity, methodsave::Bool=true,
    cache=ThreadSafeDict(), modelz=ThreadSafeDict())
    # Now let's take the formula and apply them
    formulae = OrderedDict() 
    formulae["ca1pfc"] = formula_method(df, cells, "CA1");
    formulae["pfcca1"] = formula_method(df, cells, "PFC");
    @assert !isempty(first(values(formulae)))
    glmsets = []
    for indep in ("ca1pfc", "pfcca1"), relcyc in -opt["cycles"]:opt["cycles"], 
        f in formulae[indep]
         push!(glmsets, (indep, relcyc, f))
     end
     prog = Progress(length(glmsets); desc="GLM spike counts")
     #= Threads.@threads =# for (indep, relcyc, f) in glmsets
        if unitwise
            unitstr = !isempty(unitrep) ?
                replace(string(f.lhs),unitrep...) : string(f.lhs)
            unit=parse(Int, unitstr)
            key = (;unit,indep, relcyc)
        else
            key = (;indep, relcyc)
        end
        try
           dx, dy = dx_dy[relcyc]
           cols = [string(ff) for ff in f.rhs]
           XX = Matrix(dx[!,cols])
           y  = Vector(dy[!,string(f.lhs)])
           misses = (!).(ismissing.(y))
           XX, y = xtrans.(XX[misses,:]), ytrans.(Int.(y[misses]))
           @infiltrate
           modelz[key] = if glmtool == :matlab
               if key ∈ keys(modelz) && 
                   modelz[key] isa NamedTuple
                       continue
               else
                    glm_matlab(XX, y, Dist)
               end
            elseif glmtool == :glmjl
                glm_glmjl(XX, y, Dist)
            elseif glmtool == :mlj
                glm_mlj(XX, y, Dist)
            end
           # TODO better research glmnet -- do i need to link manually -- its
           # negative output for a neural firing prediction
           # models[(;indep, relcyc, unit)] = m =  glmnetcv(XX, y, Dist)
           cache[key] = (;XX, y)
        catch exception
                @infiltrate
            modelz[(;indep, relcyc, unit)] = exception
            print("Error on $key, exception = $exception")
            sleep(0.1)
        end
        next!(prog)
     end
     @info "saving spikcount"
     if isdefined(Main, :storage); close(Main.storage); end
     storage = jldopen(fn, "a")
     try
        runname in keys(storage) ? delete!(storage, runname) : nothing
        storage[runname] = modelz
        storage[runname * "_" * "$glmtrick"] = modelz
     catch exception
        @warn exception
     else
        @info "saved"
     finally
         if isdefined(Main, :storage); close(storage); end
     end
    return modelz, cache
end

function initorget(key)
    jldopen(path_iso(opt;append="_cyclewise"), "r") do storage
        print("Possible keys", keys(storage))
        if key in keys(storage)
            storage[key]
        else
            ThreadSafeDict()
        end
    end
end

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF HAS isolated spike per cell  🔺
# ========================
model_cellhasiso, cacheiso = initorget("model_cellhasiso_mlj"), ThreadSafeDict()
run_glm!("model_cellhasiso_mlj", construct_predict_isospikecount, 
    :mlj; Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""], 
    ytrans=x->Float64(x>0), xtrans=x->Float64(x),
    modelz=model_cellhasiso, cache=cacheiso
)

model_cellhasiso_matlab = initorget("model_cellhasiso_matlab")
run_glm!("model_cellhasiso_matlab", construct_predict_isospikecount, 
    :matlab; Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""], 
    ytrans=x->Float64(x>0), xtrans=x->Float64(x),
    modelz=model_cellhasiso_matlab)

jldopen(fn, "a") do storage
    keys(storage)
end

model_cellhasiso = run_glm!("model_cellcountiso", 
    construct_predict_isospikecount, :mlj; 
    Dist=Distributions.Poisson(), 
    unitwise=true, unitrep=["_i"=>""],
    ytrans=x->Float64(x), xtrans=x->Float64(x)
)

commit_cycwise_vars()

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF SPIKE COUNTS per cell  🔺
# ========================
model_spikecount, cache = run_glm!("model_spikecount", construct_predict_spikecount,
    :mlj, Dist=Distributions.Poisson(), unitwise=true, 
    xtrans=x->Float64(x), ytrans=x->Float64(x))

commit_cycwise_vars()

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF HAS_ISO at all
# ========================
model_spikecount = run_glm("model_spikecount", construct_predict_iso,
    :mlj, Dist=Distributions.Binomial(), unitwise=false, 
    xtrans=x->Float64(x), ytrans=x->Float64(x))

# TODO add iso_count

commit_cycwise_vars()

#    _  _     ____  _       _   _   _             
#  _| || |_  |  _ \| | ___ | |_| |_(_)_ __   __ _ 
# |_  ..  _| | |_) | |/ _ \| __| __| | '_ \ / _` |
# |_      _| |  __/| | (_) | |_| |_| | | | | (_| |
#   |_||_|   |_|   |_|\___/ \__|\__|_|_| |_|\__, |
#                                           |___/ 

Plot.setfolder("isolated", "glm")
Plot.setappend((;animal=opt["animal"], day=opt["day"], tet=opt["tet"]))

grabfield(x::String) = (k,v)->if v isa Exception || v === nothing
    NaN
else
    try v[x]
    catch NaN end
end
grabfield(x::Vector{String}) = (k,v)->if v isa Exception || v === nothing
    Dict(k=>NaN for k in x)
else
    try Dict(k=>v[k] for k in x)
    catch; Dict(k=>NaN for k in x); end
end

D = to_dataframe(model_cellhasiso, grabfield("adjr2"))
Dsum = sort(combine(groupby(D, [:relcyc, :indep]), 
    :value=> (x->nanmean(collect(DIutils.skipnan(x)))) =>:value), [:indep,:relcyc])
p = begin
    @df @subset(Dsum,:indep .== "pfcca1") begin
        plot(:relcyc, :value, alpha=0.5, label="pfc->ca1")
    end
    @df @subset(Dsum,:indep .== "ca1pfc") begin
        plot!(:relcyc, :value, alpha=0.5, label="ca1->pfc")
    end
end
P=[(@df d scatter(:relcyc, :value, alpha=0.5, ylabel="r2", tickfontsize=3,
        c=:indep[1] .== "pfcca1" ? :red : :blue))
    for d in groupby(D, :unit)]
using Blink, Interact
ui = @manipulate for i in eachindex(P)
    P[i]
end; w=Window(); body!(w,ui)

D = to_dataframe(model_cellhasiso, grabfield("mae"))
Dsum = sort(combine(groupby(D, [:relcyc, :indep]), 
    :value=> (x->nanmean(collect(DIutils.skipnan(x)))) =>:value), [:indep,:relcyc])
p = begin
    @df @subset(Dsum,:indep .== "pfcca1") begin
        plot(:relcyc, :value, alpha=0.5, label="pfc->ca1")
    end
    @df @subset(Dsum,:indep .== "ca1pfc") begin
        plot!(:relcyc, :value, alpha=0.5, label="ca1->pfc")
    end
end

D = to_dataframe(model_cellhasiso, grabfield("mae"))



plot(
    plot(
        glmplot(model_isocount, "pfcca1", grabfield, label="iso count"),
        glmplot(model_hasiso,   "pfcca1", grabfield, label="iso occurance");
        link=:y,
        layout=grid(2,1),
        title= "pfc -> ca1"
    ),
    plot(
        glmplot(model_isocount, "ca1pfc", grabfield, label="iso count"),
        glmplot(model_hasiso,   "ca1pfc", grabfield, label="iso occurance");
        link=:y,
        layout=grid(2,1),
        title="ca1 -> pfc"
    ),
    ylims=(0,0.4),
    size=(600, 900)
)
Plot.save("iso, count and occur, date=$(Dates.now())")


