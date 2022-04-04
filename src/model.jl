module model
    using DrWatson
    using StatsBase
    using ProgressMeter
    using DataFrames
    using NaNStatistics
    using ..field
    import ..table
    include("utils.jl")
    include("utils/SearchSortedNearest.jl/src/SearchSortedNearest.jl")

    """
    reconstruction(jointdensity, marginalrate) -> R̂

    Example: reconstruct the response to `p(x, y, Γ)` under marginal `p(Γ)`
    """
    function reconstruction(jointdensity::AbstractArray,
            marginalrate::AbstractArray; kws...)
        dims = setdiff(1:ndims(jointdensity), findall(size(marginalrate).==1))
        #println("dims=$dims")
        A = jointdensity .* marginalrate
        B = jointdensity
        for dim ∈ dims
            A =  nansum(A; dims=dim) 
            B =  nansum(B; dims=dim)
        end
        return utils.squeeze(A./B)
    end

    """
    reconstruction_error(joint, reconstructed_joint) -> ϵ
    """
    function reconstruction_error(joint::AbstractArray, 
                                  reconstructed_joint::AbstractArray)
        Δ = diff(extrema(joint))
        ϵ = ((joint - reconstructed_joint).^2) ./ Δ
    end

    """
    data(spikes, behavior; grid, props)

    gets a form that represents the data

    get the D of spikes per behavioral Δt
    """
    function data(spikes, behavior; grid=nothing, props=nothing,
            toleranceₑ=0.0, as="group", fixbigfact=true)
        spikes = copy(spikes)

        # Map each spike to its behavioral bin
        nearest_points = SearchSortedNearest.searchsortednearest
        spikes[!,"nearest"] = nearest_points.([behavior.time], spikes.time)

        # Filter out spikes that do not match the tolerance
        Δϵ = abs.(behavior.time[spikes.nearest] .- spikes.time)
        spikes = spikes[Δϵ .< toleranceₑ, :]

        # Split by tetrode and behavioral bin ⟶   get D!
        D = combine(groupby(spikes, [:unit, :nearest], sort=true), nrow)
        D = unstack(D, :unit, :nrow)
        for (c, name) in zip(1:size(D, 2), names(D))
            col = D[!, c]
            D[ismissing.(col), c] .= 0
            if fixbigfact && name != "nearest"
                D[col .> 20, c] .= 20 # maximumal non-big() factorial
            end
        end

        # Let's lookup the linear index into a field for every behavioral
        # time
        if grid != nothing && props != nothing
            # Lookup each of the properties that we need to look into the grid fro
            properties = behavior[D[!, "nearest"], props]
            # From there, we have to find where each point in those proeprties
            # lives inside of the grid
            @assert length(grid) == length(props) "UH oh, look at the length of your passed in objects"
            properties = Matrix(properties)
            behavior[!, "fieldIndex"] = zeros(Int64, size(behavior,1))
            @showprogress for (i,row) in enumerate(eachrow(properties))
                W = fit(Histogram, Tuple(([r] for r in row)), grid).weights
                w = findall(vec(W) .== 1)
                if isempty(w) || w[1] == 0
                    continue
                end
                behavior[i, "fieldIndex"] = w[1]
            end
            D[!, "fieldIndex"] = behavior[D[!,"nearest"], "fieldIndex"]
        end

        return D
    end

    """
    probability

    grabs the probability of the models given the data at each time
    Σ probabilityₜ(dataₜ)
    """
    function probability(data::DataFrame, fieldmodels::Dict)
        out = Dict(name=>probability(data, vec(model), name) 
                    for (name,model) in fieldmodels)
        return out
    end
    function probability(data::DataFrame, fieldmodels::AbstractArray,
            units::AbstractVector)
        fieldIndex = data[!,:fieldIndex]
        probs = []
        for (model, unit) in zip(model, units)
            d = d[!, string(unit)]
            push!(probs, probability(d, model, fieldIndex))
        end
        return probs
    end
    function probability(data::DataFrame, fieldmodels::AbstractArray,
            name::NamedTuple)
        probability(data, fieldmodels, name.unit)
    end
    function probability(data::DataFrame, fieldmodel::AbstractVector,
            name::Int64)
        if string(name) in names(data)
            out = probability(data[!, string(name)], fieldmodel, data[!,:fieldIndex])
        else
            #println("Name=$name not in $(names(data))")
            out = nothing
        end
        return out
    end
    function probability(data::DataFrame, fieldmodels::AbstractArray,
            name::Int64)
        probability(data, fieldmodels[!, name], name)
    end
    function probability(counts::AbstractVector,
            fieldmodel::AbstractVector, fieldIndex::AbstractVector)
        goodInds = fieldIndex .> 0
        counts, fieldIndex = counts[goodInds], fieldIndex[goodInds]
        λ, k = fieldmodel[fieldIndex], counts
        return poisson.(k, λ)
    end
    function poisson(k::Real, λ::Real)
        out = (exp(-λ) * λ^k)/factorial(k)
        return out
    end

    """
     poisson_model : Creates a poisson_model from a field

    P⃗{x; λ}, where the model is valid ⟺   x is a count
    drawn from a single timebin Δt

    returns a function of field dimensions (to assign
    likelihood of each point)

    if grid == nothing
        then we assume users are passing in indices into our field distribution
    else
        users are passing a point, where we have to lookup the closted match
        in the grid
    end

    # modes
    K (give number of spikes), lambda is a Scalar
    lookup, lookup a behavior value of a spike, lambda is an Array
    lookup * K, K is a scalar and lambda is an Array
    """
    function poisson(k::Vector, λ::Vector, data::Tuple=nothing,
            grid::Tuple=nothing; lookup_K=true)
        if grid isa Vector && grid[1] isa Vector
            grid = Tuple(grid)
        end
        𝐊 = grid_lookup_func(grid)
        out = (exp.(-λ) .* λ.^(k * 𝐊(data)))./factorial.(k * 𝐊(data))
        return out
    end

    """
    grid_lookup_func

    creates a function that looks up field activity within a grid to
    and returns a count on that grid indicating number of events found
    per state.

    in essence, 

    it should take a point 𝐱 ∈ ℝᵈ and it will lookup it's location in
    the collection of edges 𝐆 ∈ ℝᵐ×ᵈ

    approach, we either:
    - fit() a histogram model from statsbase
    - nanstats
    """
    function grid_lookup_func(grid)
        lookup(data) = fit(Histogram, data, grid).weights
    end

    function run(spikes, beh, fieldobj, props; type="")
        if type==""
            @warn "Blank type"
        end
        # Ready the data
        D = model.data(spikes, beh; grid=fieldobj.cgrid, props=props)
        # Acquire the probabilities fields | data
        likelihood = model.probability(D, fieldobj.hist)
        # Convert to dataframe
        likelihood = field.to_dataframe(likelihood,
                                        other_labels=Dict(:type=>type),
                                        name="prob")
        likelihood[!,"logprob"] = log10.(likelihood.prob)
        table.naninf_to_missing!(likelihood, [:prob, :logprob])
        return likelihood
    end

    module plot
        using DataFrames
        using Statistics
        using Gadfly, StatsPlots, DrWatson

        function individual_poisson_mean_unitarea(likelihood; field=:prob)
            meanlike = combine(groupby(dropmissing(likelihood), [:unit, :area]), 
                               field=>mean)
            field = String(field)
            meanarealike = combine(groupby(meanlike, :area),
                                   Symbol(field*"_mean")=>mean)
            Gadfly.set_default_plot_size(20cm, 10cm)
            p = [Gadfly.plot(meanlike, color=:unit, x=:unit,
                             y=Symbol(field*"_mean"), Geom.bar),
                 Gadfly.plot(meanarealike, color=:area, x=:area,
                             y=Symbol(field*"_mean_mean"),
                             Geom.bar(position=:stack)),
                 Gadfly.plot(meanlike, color=:area, x=:area,
                             y=Symbol(field*"_mean"),
                             Geom.bar(position=:dodge))]
            hstack(p...)
        end

        function comparison_binaryop_striplot_unitarea(df; label="tmp", hline=[])
            ylabel = replace(label, "-div-"=>"/")
            p = @df df Plots.scatter(:areaNum, :prob, group=:area, ylabel=ylabel, legend=:bottomleft, alpha=0.6, xticks=((2,1),["CA1","PFC"]))
            Plots.hline!(hline; c=:black, linestyle=:dash, label="line of similarity")
            savefig(p, plotsdir("fields","poisson_model", "$label.pdf"))
            savefig(p, plotsdir("fields","poisson_model", "$label.png"))
            savefig(p, plotsdir("fields","poisson_model", "$label.svg"))
            p
        end

    end
    export plot

    module comparison
        using DataFrames
        using Statistics
        using GroupSlices
        using Distributions
        function binaryop(library::Dict, which_keys::Vector{String};
                by::Vector{String}=["unit","area"],
                op=-, measures::Vector{String}=["prob","logprob"])
            println("Which_keys = $which_keys")

                df = vcat((value for (key,value) in library if key in which_keys)...)

                splits = [by; "type"]
                gdf = groupby(df, splits)
                new = DataFrame()
                for (i,d) ∈ enumerate(gdf)
                    n = combine(dropmissing(d), :prob=>mean=>:prob, :logprob=>mean=>:logprob)
                    n = hcat(n, DataFrame(d[1, splits]))
                    append!(new, n)
                end

                out = nothing
                uType = unique(new.type)
                @assert length(uType) == 2 "not enough types"
                for (i,measure) ∈ enumerate(measures)
                    measure = [measure]
                    unused = setdiff(measures, measure)
                    u = unstack(new[!, Not(unused)], :type, measure[1])
                    if i == 1
                        out = u
                    end
                    measure=measure[1]
                    subtract = ((x,y) -> op(y,x))
                    out[!,measure]=select(u, uType => subtract => measure)[!,measure]
                end
                out = out[!, Not(uType)]
                out.areaNum = GroupSlices.groupslices(out.area) + rand(Uniform(-0.1,0.1), length(out.area))
                out
            end
    end
    export comparison

end
export model
