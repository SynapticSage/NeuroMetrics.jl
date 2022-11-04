module spiking

    import Munge
    using StatsBase
    import Utils: binning
    using DataFrames
    using ImageFiltering
    using DimensionalData
    import DimensionalData: DimArray
    using LazySets
    using ProgressMeter
    using Infiltrator
    using TensorToolbox
    using Table
    import Utils

    import Field: ReceptiveField

    export torate, rate_todataframe

    bindefault = 0.020 # 20 milliseconds
    gaussiandefault = bindefault * 3

    """
        nonlocality

    whether a data frame of behavior or spikes is NOT in a place field
    """
    function nonlocality(X::DataFrame, R::ReceptiveField; hull=1)::BitVector
        @assert :hullseg_grid ∈ R.metrics
        props = R.grid.props
        X = X[:, props]
        X = [element(Singleton(x)) for x in eachrow(X)]
        X .∉ VPolygon(R.metrics[:hullseg_grid][hull])
    end

    """
        locality

    whether a data frame of behavior or spikes is in a place field
    """
    function locality(X::DataFrame, R::ReceptiveField; hull=1)::BitVector
        @assert :hullseg_grid ∈ R.metrics
        props = R.grid.props
        X = X[:, props]
        X = [element(Singleton(x)) for x in eachrow(X)]
        X .∈ VPolygon(R.metrics[:hullseg_grid][hull])
    end


    """
        tocount

    getting spike count matrix/tensor from DF of spikes
    """
    function tocount(spikes::DataFrame, dims=:unit; binsize=bindefault, grid=nothing, kws...)
        if grid === nothing
            grid = minimum(spikes.time):binsize:maximum(spikes.time)
        end
        dims = dims isa Vector ? dims : [dims]
        T =  Munge.tensorize(spikes, dims, :time)
        prog = Progress(length(T); desc="Executing count of $dims")
        M = Array{DimArray}(undef, size(T)...)
        Threads.@threads for ind in eachindex(T)
            M[ind] = tocount(T[ind]; grid, binsize, kws...) 
            next!(prog)
        end
        neuronax = T.dims
        if grid == :dynamic
            DimArray(M, neuronax...)
        else
            timeax = M[1].dims
            M = hcat(M...)
            M = matten(M, 1, [size(M,1), size(T)...])
            DimArray(M, (timeax..., neuronax...))
        end
    end

    function tocount(times::Missing; grid, gaussian=0, binsize=nothing,
            type::Union{Nothing,Type}=nothing)::DimArray
        if grid == :dynamic
            centers = []
        else
            centers = binning.edge_to_center(collect(grid))
        end
        if gaussian > 0
            type  = type === nothing ? Float32 : type
        else
            type  = type === nothing ? UInt8 : type
        end
        DimArray(zeros(type, size(centers)), Dim{:time}(centers))
    end

    function tocount(times::AbstractArray; grid, gaussian=0, binsize=bindefault,
            type::Union{Nothing,Type}=nothing)::DimArray
        if grid === nothing
            grid = minimum(times):binsize:maximum(times)
        end
        δ = grid[2] - grid[1]
        count = fit(Histogram, vec(times), grid)
        centers = binning.edge_to_center(collect(grid))

        if gaussian > 0
            type  = type === nothing ? Float32 : type
            gaussian = gaussian * 0.1/δ
            ker = Kernel.gaussian((gaussian,))
            val = convert(Vector{type}, imfilter(count.weights, ker))
        else
            type  = type === nothing ? UInt8 : type
            val = convert(Vector{type}, count.weights)
        end
        DimArray(val, Dim{:time}(centers))
    end

    
    """
        torate (w/behavior)

    if passed with behavior, it attempts to lay out bins centered at each
    behavioral sample

    see torate(spikes::DataFrame, dims) for doc of the rest of the 
    functionalities
    """
    function torate(spikes::DataFrame, beh::DataFrame, dims=:unit; 
            binning_ratio=1, kws...)
        grid = copy(beh.time)
        δ = median(diff(beh.time)) / binning_ratio
        grid .+= δ
        grid = [[grid[1]-δ]; grid]
        grid = grid .- (1/2)δ
        # Constrain to epoch periods
        epoch_periods = Table.get_periods(beh, "epoch")
        in_period = [Table.isin.(spikes.time,
                                 epoch_period.start, epoch_period.stop) 
         for epoch_period in eachrow(epoch_periods)]
        in_period = sum(in_period)
        spikes = spikes[in_period .> 0, :]
        # Make acquire spiking structure
        torate(spikes, dims; kws..., grid)
    end

    """
        torate

    get rate matrix from a dataframe of spikes per cut of the data in 
    dims=:unit
    """
    function torate(spikes::DataFrame, dims=:unit; binsize=bindefault,
                    grid=nothing, kws...)

        if grid === nothing
            grid = minimum(spikes.time):binsize:maximum(spikes.time)
        end
        dims = dims isa Vector ? dims : [dims]
        T =  Munge.tensorize(spikes, dims, :time)
        M = Array{DimArray}(undef, size(T)...)
        prog = Progress(length(T); desc="Executing count of $dims")
        for ind in eachindex(T)
            M[ind] = torate(T[ind]; grid, binsize, kws...) 
            next!(prog)
        end
        neuronax = T.dims
        if grid == :dynamic
            DimArray(M, neuronax...)
        else
            timeax = M[1].dims
            M = hcat(M...)
            M = matten(M, 1, [size(M,1), size(T)...])
            DimArray(M, (timeax..., neuronax...))
        end
    end

    torate(times::Missing; kws...)::DimArray = tocount(times; kws...)


    function torate_windowdia(times::AbstractArray; grid, windowsize::Real,
            gaussian::Real=0)
        torate_windowrad(times; grid, radius=windowsize/2, gaussian)
    end
    #function torate(times::AbstractArray; grid, windowsizes::Tuple{<:Real,<:Real})
    #end
    function torate_windowrad(times::AbstractArray; grid, radius::Real, 
                           gaussian::Real=0)
        count = binning.inside(times, grid, radius)
        ker = Kernel.gaussian((gaussian,))
        δ  = grid[2] - grid[1]
        centers = collect(grid)
        DimArray(imfilter(count ./ δ, ker),
                  Dim{:time}(centers))
    end

    function torate(times::AbstractArray; grid, gaussian=gaussiandefault,
            binsize=grid[2]-grid[1])::DimArray
        count = fit(Histogram, vec(times), grid)
        gaussian = gaussian * 0.1/binsize
        ker = Kernel.gaussian((gaussian,))
        centers = binning.edge_to_center(collect(grid))
        DimArray(imfilter(count.weights ./ binsize, ker),
                  Dim{:time}(centers))
    end

    # CELL COFIRING
    function xcorr(rate::DimArray, cell1::T , cell2::T; lags=-20:20) where 
        T<:Union{Int16,Int32,Int64}
        x, y = rate[Dim{:unit}(cell1)],
               rate[Dim{:unit}(cell2)]
        StatsBase.crosscor(x, y, lags)
    end

    function xcorr(spikes::DataFrame; lags=-200:200, kws...)
        units1 = unique(spikes.unit)
        units2 = unique(spikes.unit)
        R = torate(spikes; kws...)
        results = []
        for (cell1,cell2) in Iterators.product(units1,units2)
            push!(results, xcorr(R, cell1, cell2; lags))
        end
    end

    function nextandprev!(spikes::DataFrame)
        combine(groupby(spikes, :unit), nextandprev!)
    end
    function nextandprev!(spikes::SubDataFrame)
        sort!(spikes, :time)
        for field in [:prevt, :nextt, :prevd, :nextd, :neard]
            if field ∉ propertynames(spikes)
                spikes[!,field] = fill(NaN, size(spikes,1))
            end
        end
        for (r,row) in enumerate(eachrow(spikes))
            p, n = max(1,r-1), min(r+1,size(spikes,1))
            row.prevt = spikes[p,:time]
            row.nextt = spikes[n,:time]
            row.prevd = row.time - spikes[p,:time]
            row.nextd = spikes[n,:time] - row.time
            row.neard = min(row.prevd, row.nextd)
        end
        spikes
    end

    """
        isolated

    find isolated spikes in the manneer of Jai/Frank 2021
    """
    function isolated(spikes::DataFrame,  theta::Union{DataFrame,Nothing}; 
                      cycle=:cycle, kws...)
        if !hasproperty(spikes, Symbol(cycle))
            Utils.filtreg.register(theta, spikes; on="time", transfer=[String(cycle)])
        end
        prog = Progress(length(unique(spikes.unit)); desc="Adding isolation stats")
        combine(groupby(spikes, :unit), 
                     x->(i=isolated(x; kws...);next!(prog);i)
                     )
    end
    function isolated(spikes::SubDataFrame; N=3, thresh=8, cycle_prop=:cycle, include_samples::Bool=false)
        explore = setdiff(-N:N,0)
        if !hasproperty(spikes, :isolated)
            spikes[!,:isolated] = Vector{Union{Missing,Bool}}(missing, size(spikes,1))
        end
        if !hasproperty(spikes, :nearestcyc)
            spikes[!,:nearestcyc] = Vector{Union{Missing,Int32}}(missing, size(spikes,1))
        end
        if !hasproperty(spikes, :meancyc)
            spikes[!,:meancyc] = Vector{Union{Missing,Float32}}(missing, size(spikes,1))
        end
        if include_samples && !hasproperty(spikes, :isosamples)
            spikes[!,:meancyc] = Vector{Union{Missing,Vector}}(missing, size(spikes,1))
        end
        cycles = groupby(spikes, cycle_prop)
        #if length(cycles) > 1
            #@warn  "You only have 1 cycle"
        #end
        #(c,cycle) =  first(enumerate(cycles))
        for (c,cycle) in enumerate(cycles)
            # find the N closest cycles
            explore_cycles = unique(max.(min.(c .+ explore, [size(cycles,1)]),[1]))
            if length(explore_cycles) < length(explore)
                cycle.isolated .= false

                center_times = [abs(mean(c.time)-mean(cycle.time)) for c in cycles[explore_cycles]]
                order = sortperm(center_times)
                nearest = explore_cycles[order] # TODO
                cycle_prox = [abs(cycle[1,cycle_prop] - other_cyc[1,cycle_prop])
                              for other_cyc in cycles[nearest]]
                nearestcyc = minimum(cycle_prox)
                meancyc    = mean(cycle_prox)
                cycle.meancyc    .= meancyc
                cycle.nearestcyc .= nearestcyc

            else
                center_times = [abs(mean(c.time)-mean(cycle.time)) for c in cycles[explore_cycles]]
                order = sortperm(center_times)
                nearest = explore_cycles[order[1:N]] # TODO make this match TODO above
                cycle_prox = [abs(cycle[1,cycle_prop] - other_cyc[1,cycle_prop])
                              for other_cyc in cycles[nearest]]
                nearestcyc = minimum(cycle_prox)
                meancyc    = mean(cycle_prox)
                isolated   = meancyc > thresh
                cycle.isolated   .= isolated
                cycle.meancyc    .= meancyc
                cycle.nearestcyc .= nearestcyc
                if include_samples
                    cycle.isosamples .= [cycle_prox]
                end
            end
        end
        spikes
    end


    """
        rate_todataframe

    converts a rate dimarray (or count dimarray) to a dataframe, and populates columns labeling the
    times from other dataframes
    """
    function rate_todataframe(X::DimArray, registrant::Tuple{DataFrame, String, Vector{String}}...)::DataFrame
        Table.from_dimarrary(X, registrant...)
    end

end



