module kerneldens

    using ..Field
    using DataFrames
    using StatsBase, KernelDensity
    using KernelDensitySJ
    #using Infiltrator
    #using DrWatson
    import Utils

    function KDE(data, props; bandwidth=:silverman)
        data = dropmissing(data[:,props])
        #for col in props
        #    inds = data[:, col] .== NaN
        #    data[inds, col] .= 0
        #end
        if isempty(data) || data isa Nothing
            return nothing
        else
            fail=false
            bandwidths=nothing
            if bandwidth == :sj
                try
                    bandwidths = tuple([bwsj(data[:,prop]) for prop in props]...)
                catch e
                    fail=true
                end
                if fail || !all(bandwidths .> 0) 
                    return nothing
                end
            elseif bandwidth == :silverman
                bandwidths = nothing
            else
                bandwidths = bandwidth
            end
            if bandwidths == nothing
                return kde(tuple([convert(Vector{Float32}, x) 
                                  for x in eachcol(data[:, props])]...
                                )
                          )
            else
                return kde(tuple([convert(Vector{Float32}, x) 
                                  for x in eachcol(data[:, props])]...
                                ); bandwidth=bandwidths
                          )
            end
        end
    end


    function gridpdf(kde, grid)
        grid = [grid[prop] for prop in keys(grid)]
        pdf(kde, grid...)
    end

    """
        norm_kde_by_histcount

    normalizes kernel density estimates to have the same number of event
    counts as the histogram, so that it will accurately represent peak
    firing/ripple rates
    """
    function norm_kde_by_histcount(kde::Union{AbstractArray, Dict}, 
            hist::Union{AbstractArray,Dict}, skip_keyerror=false)
        #if typeof(hist) ≠ typeof(kde)
        #    throw(ArgumentError("type of hist should match type of kde\n...type(hist)=$(supertype(typeof(hist))) != $(supertype(typeof(kde)))"))
        #end
        if kde isa AbstractArray
            area_hist = sum(Utils.skipnan(hist))
            area_kde  = sum(Utils.skipnan(kde))
            kde = (area_hist/area_kde) .* kde;
        elseif kde isa Dict
            @inbounds for key in keys(kde)
                try
                    kde[key] = norm_kde_by_histcount(kde[key], hist[key])
                catch KeyError
                    if !skip_keyerror
                        throw(KeyError("Key=$key missing from hist"))
                    end
                end
            end
        else
            throw(ArgumentError)
        end
        return kde
    end

    function fields(data::DataFrame, beh::DataFrame;
            splitby::Union{Nothing,Vector{String},Vector{Symbol},String}="unit",
            resolution::Union{Int,Vector{Int}}=200,
            savemem::Bool=true,
            props::Vector{String}=["x","y"])

        # Grid settings
        grid = getSettings(beh, props, 
                                 resolution=resolution, settingType="kde")
        #@infiltrate
        grid_hist = getSettings(beh, props, resolution=resolution,
                                settingType="hist")
        # Behavioral distribution
        behDist = KDE(beh, props);
        behDist_hist = hist.h2d(beh, props, grid=grid_hist).weights;

        # Data distribution applying splits if applicable
        ith_group(i) = DataFrame(groups[i])
        field_of_group(i) = KDE(ith_group(i), props)
        if splitby isa Vector{Symbol}
            splitby = String.(splitby)
        end
        if all(in.(splitby, [names(data)]))
            groups = groupby(data, splitby)
            ugroups = collect(zip(sort(collect(groups.keymap),
                                       by=x->x[2])...))[1]
            ugroups = 
            [(;zip(groups.cols,ugroup)...) for ugroup in ugroups]
            spikeDist = Dict(ugroups[i] => try field_of_group(i) catch nothing end
                             for i ∈ 1:length(groups)) # TODO replace nothing with array of nan
        else
            if splitby != nothing
                @warn "splitby not empty, yet not all columns in df"
            end
            spikeDist = KDE(data[!,:], props)
        end

        # Occupancy
        behDist = pdf(behDist, grid...);
        zero_fraction = 0.0000000001
        behDist[behDist .<= zero_fraction] .= 1
        behzeroinds = behDist_hist .<= 0
        if savemem
            behDist = Float32.(behDist)
        end

        # Spike count
        if spikeDist isa Dict
            dist = Dict{typeof(ugroups[1]),Any}();
            @inbounds for i ∈ keys(spikeDist)
                if spikeDist[i] == nothing
                    dist[i] = nothing
                else
                    if savemem
                        dist[i] = Float32.(pdf(spikeDist[i], grid...))
                    else
                        dist[i] = pdf(spikeDist[i], grid...)
                    end
                end
            end
        else
            @assert dist != nothing
            dist = pdf(spikeDist, grid...)
        end

        return (kde=dist, grid=grid, occ=behDist, occzeroinds=behzeroinds)
    end

end
