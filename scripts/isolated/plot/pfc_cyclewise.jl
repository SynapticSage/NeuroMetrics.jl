# IMPORTS AND DATA
# --------
checkpoint = false
include("../imports_isolated.jl")
include("./imports_isolated.jl")
if checkpoint
    include(scriptsdir("isolated","load_cyclewise_checkpoint.jl"))
    @time include("../load_cyclewise_checkpoint.jl")
else
    @time include(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")
end
include(scriptsdir("isolated","imports_isolated.jl"))


# Filtration
pyramidal = cells[cells.meanrate .< 5, :unit]
spikecount = combine(groupby(spikes, :unit), 
    nrow=>:spikecount, :isolated=>(x->sum(collect(skipmissing(x))))=>:isospikecount)
replace!(spikecount.isospikecount, missing=>0)
replace!(spikecount.spikecount, missing=>0)
spikecount_criteria = 
    spikecount[spikecount.isospikecount .> 5 .&& spikecount.spikecount .> 50,:unit]
neurons_of_interest = intersect(pyramidal, spikecount_criteria)
spikes_sub = @subset(spikes, :unit .∈ (neurons_of_interest,))
Rdf_sub    = @subset(Rdf, :unit .∈ (neurons_of_interest,))
cells_sub  = @subset(cells, :unit .∈ (neurons_of_interest,))
@assert size(spikes,1) != size(spikes_sub,1)
spikes_sub = dropmissing(spikes_sub, :isolated)

commit_vars()
commit_cycwise_vars()

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

    # -----------------------------------------------------
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

    # -----------------------------------------------------
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


    # --------------------------------------------------------------------
    # Figure out theta cycles with theta power at least 1/2 std above mean
    # and those which are accounted for in the grid occupancy
    # --------------------------------------------------------------------
    DIutils.filtreg.register(beh, lfp, on="time", transfer=["speed"])
    upper = quantile(@subset(cycles, :speed .> 2, :amp_mean .< 300).amp_mean, 0.995)
    lower = quantile(@subset(cycles, :speed .> 2, :amp_mean .< 300).amp_mean, 0.005)
    DataFrames.transform!(cycles, [:amp_mean] => (a-> a .> lower .&& 
        a .< upper) => :theta_good, renamecols=false)
    cycles.speed = Vector{Union{Missing,Float32}}(missing, size(cycles,1))
    for (i,(start, stop)) in enumerate(zip(cycles.start, cycles.stop))
        inds = DIutils.in_range(beh.time, [start, stop])
        cycles[i, :speed] = mean(beh[inds,:speed])
    end
    # Label cyccles as having occupancy
    cycles.hasocc = (!).(ismissing.(occ.datainds))
    # Register that across dataframes
    kws = (;on="cycle", transfer=["hasocc", "theta_good", "speed"])
    DIutils.filtreg.register(cycles, Rdf;     kws...)
    DIutils.filtreg.register(cycles, Rdf_sub; kws...)


    # -----------------------------------------------------
    # Prepare to select good cycles
    # -----------------------------------------------------
    # Determine isolated cycles to look at
    iso_cycles = unique(@subset(Rdf_sub, 
                        :isolated_sum .> 0, 
                        :hasocc .== true,
                        :theta_good .== true,
                        :speed .> 2
    ).cycle)
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

# DIutils.filtreg.register(cycles, df, on="cycle", transfer=["epoch","time","bins"]);
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
histogram(df.cycs, label="df", title="Check for uniformity in cyc production")
histogram2d(df.relcycs, df.bins)
checkbins(:df);

# Remove missing vals from isolated spike count columns
for col in eachcol(df[!,[x for x in names(df) if occursin("_i", x)]] )
    replace!(col, missing=>0)
    col = convert(Vector{Union{Missing,Float64}}, col)
end

# Checkpoint
commit_cycwise_vars()
commit_cycwise_vars("df")
opt["commit"] = true
opt["overwrite"] = true

# Compare df to loaded df
    dd = initorget("df")
    sort!(df, [:relcycs, :cycs])
    sort!(dd, [:relcycs, :cycs])
    size(df), size(dd)
    df[:, cellcols(df)] .== dd[:, cellcols(dd)]
    plot(
    heatmap(log10.(Matrix(df[:, cellcols(df)])), title="df"),
    heatmap(log10.(Matrix(dd[:, cellcols(dd)])), title="dd")
    )
    histogram(df.cycs,  label="df")
    histogram!(dd.cycs, label="dd")
    dd.cycle = dd.cycs
    DIutils.filtreg.register(cycles, dd, on="cycle", transfer=["time"])
    DIutils.filtreg.register(cycles, df, on="cycle", transfer=["time"])
    histogram(df.time, label="df")
    histogram!(dd.time, label="dd")


# import CSV
# CSV.write(expanduser("~/df_sub.csv"), df_sub)
