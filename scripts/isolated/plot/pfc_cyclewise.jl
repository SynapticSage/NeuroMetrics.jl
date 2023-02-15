if !(:lfp in names(Main))
    # @time include(scriptsdir("isolated","load_isolated.jl"))
    include("../load_isolated.jl")
end

# Get dataframe of R
val = :value
Rdf = DataFrame(R; name=val)
# Set cycle via start cycle times
cycles.time = cycles.start
DIutils.filtreg.register(cycles, Rdf, on="time", 
                         transfer=["cycle"], 
                         match=:prev)


# Figure out cycles with some number of isolated spikes
cycstat = combine(groupby(@subset(spikes, :area .== "CA1"), 
                :cycle), 
        :isolated => sum, 
        :unit => (n->unique(n)) => :diversity,
        [:isolated, :unit] => ((i,u)->(length(unique(u[i.==true])))) => :isodiversity,
        [:isolated, :unit] => ((i,u)->(length(unique(u[i.==false])))) => :adjdiversity,
       )
dropmissing!(cycstat, :cycle)


fn = path_iso(opt; append="_cyclewise")
if !isfile(fn)

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
    Rdf_cycles = groupby(Rdf, [:cycle])
    Rdf_isocycles = unique(@subset(Rdf, :isolated_sum .> 0).cycle)
    uArea = unique(cells.area)
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

    df, cyc_error = Vector{Union{Missing,DataFrame}}(missing, length(Rdf_isocycles)), 
                    Dict()
    Infiltrator.clear_disabled!()
    prog = Progress(length(Rdf_isocycles); desc="cycle df")
    Threads.@threads for (i,cyc) in collect(enumerate(Rdf_isocycles))
        try
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
            df[i] = hcat(DataFrame([cycs,relcycs],[:cycs,:relcycs]), 
                         vcat(U...; cols=:union))

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

# PREPROCESS
begin
    # Helper functions :hand:
    using GLM
    function construct_predict_spikecount(df, independent_area="CA1";
            other_vars=[], other_ind_vars=[])
        uArea = unique(cells.area)
        @assert length(uArea) == 2 "Only supports two area dataframes"
        dependent_area = setdiff(uArea, [independent_area])

        dep_neurons = @subset(cells,:area .==dependent_area).unit
        ind_neurons = @subset(cells,:area .==independent_area).unit
        filter!(n->string(n) ∈ names(df), dep_neurons)
        filter!(n->string(n) ∈ names(df), ind_neurons)
        
        formulae = []
        for nd in dep_neurons
            ni = first(ind_neurons)
            independents = GLM.Term(Symbol(ni))
            for ni in ind_neurons[2:end]
                independents += GLM.Term(Symbol(ni)) 
            end
            formula = GLM.Term(Symbol(nd)) ~ independents
            push!(formulae, formula)
        end
        formulae
    end
    function construct_predict_iso(df, independent_area="CA1", type=:has;
            other_vars=[], other_ind_vars=[])
        uArea = unique(cells.area)
        @assert length(uArea) == 2 "Only supports two area dataframes"
        ind_neurons = @subset(cells,:area .==independent_area).unit
        filter!(n->string(n) ∈ names(df), ind_neurons)
        formulae = []
        ni = first(ind_neurons)
        independents = GLM.Term(Symbol(ni))
        for ni in ind_neurons[2:end]
            independents += GLM.Term(Symbol(ni)) 
        end
        if type == :has
            formula = GLM.Term(:has_iso) ~ independents
        elseif type == :count
            formula = GLM.Term(:isolated_sum) ~ independents
        end
        push!(formulae, formula)
        formulae
    end

end

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
# ========================
#
# Requirements
# ------------
# We next identified another control set of theta cycles when the CA1 cell did
# not spike. These control cycles were matched for animal speed, movement
# direction and location (see “Cycle matching”). We created a model for each
# CA1 cell to determine whether PFC spiking activity can distinguish between
# cycles with or without isolated spiking in a time window relative to the
# cycle with isolated spiking. We first modeled using activity in the 12 cycles
# previous to the cycle with isolated activity and then the 12 cycles after the
# cycle with isolated activity. We ensured that no other isolated activity
# occurred in this window was used for prediction. A 4-cycle bin size was used
# for grouping PFC activity since PFC activity shows relatively long
# autocorrelation times.

# Now let's take the formula and apply them
using Lasso, Distributions, ThreadSafeDicts
formulae, models   = OrderedDict(), ThreadSafeDict()
formulae["ca1pfc"] = construct_predict_spikecount(df, "CA1");
formulae["pfcca1"] = construct_predict_spikecount(df, "PFC");
@assert !isempty(first(values(formulae)))
glmsets = []
for indep in ("ca1pfc", "pfcca1"), relcyc in -8:8, f in formulae[indep]
    push!(glmsets, (indep, relcyc, f))
end
using ThreadSafeDicts
prog = Progress(length(glmsets); desc="GLM spike counts")
Threads.@threads for (indep, relcyc, f) in glmsets
    try
    y, X = modelcols(apply_schema(f, schema(f, df)), df) 
    PoisMod = fit_mle(Poisson, y)
    models[(;indep, relcyc)] = fit(GLM.GeneralizedLinearModel, X, min.(y,[1]), PoisMod)
    catch
        models[(;indep, relcyc)] = nothing
    end
    next!(prog)
end

pvals = OrderedDict(k=>ftest(v) for (k,v) in models)

