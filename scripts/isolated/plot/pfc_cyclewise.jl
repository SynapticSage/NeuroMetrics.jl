if !(:lfp in names(Main))
    # @time include(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")
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

# Checkpoint
fn=replace(path_iso(opt["animal"], opt["day"]),
        ".jld2"=>"_cyclewise.jld2")
rm(fn)
jldsave(fn; df)

# Group the dataframe into subdataframes
# by relative cycle and area
dfg  = groupby(df, [:relcycs])
#dfgc = groupby(df, [:relcycs, :area, :unit])

# Helper functions :hand:
using GLM
function getarea(df::DataFrame, area::String)
    df = @subset(df, :area .== area)
    nonmissing_cols = [!all(ismissing.(col)) for col in eachcol(df)]
    df[!,nonmissing_cols]
end
Base.adjoint(s::String) = s
function construct_area_formula(df, independent_area="CA1";
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

# ========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
# ========================
# Now let's take the formula and apply them
formulae, models = OrderedDict(), OrderedDict()
formulae["ca1pfc"] = construct_area_formula(df, "CA1");
formulae["pfcca1"] = construct_area_formula(df, "PFC");
for indep in ("ca1pfc", "pfcca1"), relcyc in -8:8, f in formulae[indep]
    models[indep] = lm(f, @subset(df, :relcycs .== relcyc))
end


