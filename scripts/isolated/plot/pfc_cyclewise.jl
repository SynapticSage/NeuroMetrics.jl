if !(:lfp in names(Main))
    # @time Uinclude(scriptsdir("isolated","load_isolated.jl"))
    @time include("../load_isolated.jl")
end

import DI, DIutils.Table

# Get dataframe of R
Rdf = DataFrame(R; value)
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

# Annotate spikes and Rdf
dropmissing!(Rdf, :isolated_sum)
Rdf_cycles = groupby(Rdf, :cycle)
Rdf_isocycles = unique(@subset(Rdf, :isolated_sum .> 0).cycle)

