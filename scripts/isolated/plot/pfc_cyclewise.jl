if !(:lfp in names(Main))
    # @time include(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")
end

import DI, DIutils.Table
val = :value

# Get dataframe of R
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
Rdf_cycles = groupby(Rdf, :cycle)
Rdf_isocycles = unique(@subset(Rdf, :isolated_sum .> 0).cycle)

kws=(;label="")
plot(
     (@df cycstat histogram(:isolated_sum ;xlabel="iso spikes emitted", kws...)),
     (@df cycstat histogram(:isodiversity;xlabel="# of uniq iso cells", kws...)),
     (@df cycstat histogram(:adjdiversity;ylabel="# of uniq adj cells", kws...)),
     Plot.blank((plot();Plots.annotate!(0.5,0.5,text("Cycle statistics",14))), visible=false, size=(100,50)),
     layout=grid(2,2))

# Get the firing rates per cycle
indexers = [:time,:isolated_sum, :area]
# Dcyc = Dict(
#             :fr  => [],
#             :cyc => [],
#             :relcyc => [],
#             (x=>[] for x in indexers)...
#            )
df, cyc_error = [], Dict()
begin

    using Infiltrator
    prog = Progress(length(Rdf_isocycles); desc="cycle df")
    #= Threads.@threads =# for cyc in Rdf_isocycles
        
        try
             # Address cycles of interest
             🔑s = [(;cycle=cyc) for cyc in UnitRange(cyc-8, cyc+8)]

            # Grab each cycle of activity
            U = [begin
                     @infiltrate
                 u = unstack(Rdf_cycles[🔑], indexers, :unit, val)
                 combine(u, Not(:time) .=> [mean], renamecols=false)
             end
                for 🔑 in 🔑s if 🔑 in keys(Rdf_cycles)]
            
            cycs = hcat(([[🔑.cycle, relcyc] for (🔑,relcyc) in zip(🔑s, UnitRange(-8,8)) 
                          if 🔑 in keys(Rdf_cycles)])...)'

            # Added df to list
            push!(df, hcat(DataFrame(cycs,[:cycs,:relcycs]), vcat(U...)))
            next!(prog)
            # push!(Dcyc[:fr], U[:, length(indexers)+1:end])
            # push!(Dcyc[:cyc], cycs)
            # push!(Dcyc[:relcyc], relcycs)
            # for indexer in indexers
            #     push!(Dcyc[indexer], U[:, indexer])
            # end
        catch exception
            cyc_error[cyc] = exception
        end

    end
    @info cyc_error
    df = vcat(df...)
end

# Checkpoint
fn=replace(path_iso(opt["animal"], opt["day"]),
        ".jld2"=>"_cyclewise.jld2")
fn = open(fn,"r")
jldsave(fn; compress=true, df, cyc_error)

# Group the dataframe into subdataframes
# by relative cycle and area
dfg = groupby(df,  [:relcycs, :area])
dfgc = groupby(df, [:relcycs, :area, :unit])

# Let

using GLM
function construct_area_formula(df, area)
end


