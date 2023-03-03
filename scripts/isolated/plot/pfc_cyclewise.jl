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
Rdf_cycles  = groupby(Rdf, [:cycle])
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

# 
#    _  _     ____                                
#  _| || |_  |  _ \ _ __ ___ _ __   __ _ _ __ ___ 
# |_  ..  _| | |_) | '__/ _ \ '_ \ / _` | '__/ _ \
# |_      _| |  __/| | |  __/ |_) | (_| | | |  __/
#   |_||_|   |_|   |_|  \___| .__/ \__,_|_|  \___|
#                           |_|                   
#                  _       _       _   _          _        
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
iso_cycles = unique(@subset(Rdf, :isolated_sum .> 0, 
                                          :hasocc .== true).cycle).cycles
match_cycles!(cycles, Rdf, occ; matches=opt["matches"], iso_cycles)

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
# (each iso/noniso cycle plus precedents and antecedents)
df = df_FRpercycle_and_matched(cycles, Rdf_cycles, beh; 
                               indexers, cycrange=opt["cycles"])

# Checkpoint
commit_cycwise_vars()

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
using Distributed
"""
Shortcut function that handles applying my statsmodels formulae to each
cut of the data andn running GLM with the chosen `glmtool`
"""

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF HAS isolated spike per cell  🔺
# ========================
model_cellhasiso, cacheiso, shuffle_cellhasiso = 
                             initorget("model_cellhasiso"), 
                             ThreadSafeDict(),
                        initorget("model_cellhasiso_mlj_shuffle";obj=Dict())
pos = (df, "model_cellhasiso_mlj", construct_predict_isospikecount, :mlj)
kws = (Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""],
    ytrans=x->Float64(x>0), xtrans=x->Float64(x), modelz=model_cellhasiso,
    cache=cacheiso)
run_glm!(pos...;kws...)
run_shuffle!(pos...;kws..., shuffle_models=shuffle_cellhasiso, shufcount=20)
commit_cycwise_vars()

model_cellhasiso_matlab = initorget("model_cellhasiso_matlab")
run_glm!(df, "model_cellhasiso_matlab", construct_predict_isospikecount, 
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
model_spikecount = run_glm!(df, "model_spikecount", construct_predict_iso,
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


modelsort(x) = sort(OrderedDict(x), by= q->[q.unit, q.relcyc])
It = Iterators.product(zip(modelsort.((model_cellhasiso, model_spikecount)), 
            ("cell has iso", "spike count")), ["adjr2", "mae"])

P = Dict()
for ((model, name), measure) in It
    D = to_dataframe(model, grabfield(measure))
    Dsum = sort(combine(groupby(D, [:relcyc, :indep]), 
        :value=> (x->nanmean(collect(DIutils.skipnan(x)))) =>:value),
            [:indep,:relcyc])
    title = "$name $measure"
    p = begin
        @df @subset(Dsum,:indep .== "pfcca1") begin
            plot(:relcyc, :value, alpha=0.5, label="pfc->ca1")
        end
        @df @subset(Dsum,:indep .== "ca1pfc") begin
            plot!(:relcyc, :value, alpha=0.5, label="ca1->pfc", ylabel=measure)
        end
    end
    push!(P, title=>p)
end
plot(values(P)...)

using Blink, Interact, Interpolations
function interp(D)
    D = copy(D)
    x,y,z = 1:size(D,1), D.value, D.relcyc
    notmissing = (!).(isnan.(y))
    xn, yn, zn = x[notmissing], y[notmissing], z[notmissing]
    itp = linear_interpolation(xn, yn)
    itpr = linear_interpolation(xn, zn)
    x = unique(clamp.(x, xn[1], xn[end]))
    D[x,:value] = (itp[x,:])
    D[x,:relcyc] = (itpr[x,:])
    D = D[(!).(isnan.(D.value)), :]
    D
end

for ((model, nm), measure) in It

    tle="$nm $measure"
    D = sort(to_dataframe(model, 
                 grabfield(measure)), [:relcyc])
    P=[(#d=interp(d);
        @df d (
        scatter(:relcyc, :value; alpha=0.5, ylabel=measure, 
        title=tle, linestyle=:solid, tickfontsize=6, label=replace(:indep[1],
            "pfcca1"=>"pfc ⟶ ca1",
            "ca1pfc"=>"ca1 ⟶ pfc"),
    c=:indep[1] .== "pfcca1" ? :red : :blue)); 
    #@df d plot!(:relcyc, :value)
    )
    for d in groupby(D, :unit)]
    ui = @manipulate for i in eachindex(P)
        P[i]
    end
    w=Window(); body!(w,ui)

end


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

