module spiking

    using Base: test_success
    using StatsBase
    using DataFrames
    using ImageFiltering
    using DimensionalData
    import DimensionalData: DimArray
    using LazySets
    using ProgressMeter
    using Infiltrator
    using TensorToolbox

    import DIutils: binning
    using DIutils: Table
    import DIutils
    import ..Munge
    import ...Field: ReceptiveField

    export torate, rate_todataframe

    bindefault = 0.020 # 20 milliseconds
    gaussiandefault = bindefault * 3
    mag_constant = 0.5

    SymStr          = Union{Symbol, String, Int} # allowed column names

    """
        nonlocality(X::DataFrame, R::ReceptiveField; hull=1)::BitVector

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
        locality(X::DataFrame, R::ReceptiveField; hull=1)::BitVector

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
        tocount (w/behavior)

    if passed with behavior, it attempts to lay out bins centered at each
    behavioral sample

    see torate(spikes::DataFrame, dims) for doc of the rest of the 
    functionalities

    # Arguments
    - `spikes::DataFrame`: dataframe of spikes
    - `beh::DataFrame`: dataframe of behavior
    - `dims=:unit`: dimensions to tensorize over
    - `binning_ratio=1`: ratio of binning to behavioral sampling
    - `downsample=nothing`: downsample behavioral data by this factor
    - `kws...`: keyword arguments to pass to torate(spikes::DataFrame, dims)
    """
    function tocount(spikes::DataFrame, beh::DataFrame, 
    dims::Union{T,Vector{T}} where T <: SymStr=:unit; 
    downsample::Union{Int,Nothing}=nothing,
            binning_ratio=1, kws...)
        beh = downsample === nothing ? beh : beh[1:downsample:end, :]
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
        tocount(spikes, dims; kws..., grid)
    end


    """
        tocount

    getting spike count matrix/tensor from DF of spikes
    """
    function tocount(spikes::DataFrame, dims=:unit; binsize=bindefault, grid=nothing, kws...)
        if grid === nothing
            grid = minimum(spikes.time):binsize:maximum(spikes.time)
        end
        @assert !(dims isa DataFrame)
        dims = dims isa Vector ? dims : [dims]
        T =  Munge.tensor.tensorize(spikes, dims, :time)
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
            gaussian = gaussian * mag_constant/δ
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

    # Arguments
    - `spikes::DataFrame`: dataframe of spikes
    - `beh::DataFrame`: dataframe of behavior
    - `dims=:unit`: dimensions to tensorize over
    - `binning_ratio=1`: ratio of binning to behavioral sampling
    - `downsample=nothing`: downsample behavioral data by this factor
    - `kws...`: keyword arguments to pass to torate(spikes::DataFrame, dims)
    """
    function torate(spikes::DataFrame, beh::DataFrame, 
    dims::Union{T,Vector{T}} where T <: SymStr=:unit; 
    downsample::Union{Int,Nothing}=nothing,
            binning_ratio=1, kws...)
        beh = downsample === nothing ? beh : beh[1:downsample:end, :]
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
        torate(spikes::DataFrame, dims=:unit; binsize=bindefault, grid=nothing,
                kws...)

    get rate matrix from a dataframe of spikes per cut of the data in 
    dims=:unit
    """
    function torate(spikes::DataFrame, 
    dims::Union{T, Vector{T}} where T<:SymStr=:unit; binsize=bindefault,
    grid=nothing, kws...) 

        if grid === nothing
            grid = minimum(spikes.time):binsize:maximum(spikes.time)
        end
        dims = dims isa Vector ? dims : [dims]
        T =  Munge.tensor.tensorize(spikes, dims, :time)
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


    """
        torate_windowdia(times::AbstractArray; grid, windowsize::Real,
            gaussian::Real=0)
    """
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

    """
        torate(times::AbstractArray; grid, gaussian=0, binsize=bindefault)

    get rate matrix from a vector of times

    # Arguments
    - `times::AbstractArray`: vector of times
    - `grid`: grid to bin times on
    - `gaussian`: gaussian smoothing parameter
    - `binsize`: binsize to use if grid is not provided

    # Returns
    - `DimArray`: rate matrix
    """
    function torate(times::AbstractArray; grid, gaussian=gaussiandefault,
            binsize=grid[2]-grid[1])::DimArray
        count = fit(Histogram, vec(times), grid)
        gaussian = gaussian * mag_constant/binsize
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

    """
        xcorr(spikes::DataFrame; lags=-200:200, kws...)

    get cross correlation matrix from a dataframe of spikes per cut of the data
    in dims=:unit

    # Arguments
    - `spikes::DataFrame`: dataframe of spikes
    - `lags=-200:200`: lags to compute xcorr over
    - `kws...`: keyword arguments to pass to torate(spikes::DataFrame, dims)
    """
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

    update
    """
    function prepiso(spikes::AbstractDataFrame,  
                      theta::Union{AbstractDataFrame,Nothing}; 
                      cycle=:cycle, refreshcyc=true, cells=nothing,
                      beh=nothing, ripples=nothing, immobility_thresh=2,
                      matchtetrode::Bool=false, kws...)

        if spikes.animal|>unique|>length > 1
            @error "Can only run isolated on one animal at a time"
        end

        if refreshcyc || Symbol(cycle) ∉ propertynames(spikes)
            spikes[!,String(cycle)] = 
                Vector{Union{Float32,Missing}}(missing, size(spikes,1))
        elseif refreshcyc
            spikes[!,String(cycle)] .= missing
        end

        added_lfp_field = true
        if String(cycle) ∉ propertynames(theta)
            theta[!,String(cycle)] = theta[!,:cycle]
        end

        if immobility_thresh > 0
            if beh !== nothing
                DIutils.filtreg.register(beh, spikes; on="time", 
                                                transfer=["speedsmooth"])
            end
            if :smoothvel ∈ propertynames(spikes)
                println("Excluding immobility")
                spikes = subset(spikes, :speedsmooth => v->v.<immobility_thresh, 
                view=true, skipmissing=true)
            end  
        end

        if ripples !== nothing
            println("Excluding ripples")
            spikes_excluded = DIutils.in_range(spikes.time, ripples)
            spikes = spikes[.!spikes_excluded,:]
        end

        if !matchtetrode
            DIutils.filtreg.register(theta, spikes; on="time", 
                                     transfer=[String(cycle)])
        else
            @assert cells !== nothing
            cell_to_tet = Dict(cell=>tet for (cell,tet) in 
                                zip(cells.unit, cells.tetrode))
            G  = groupby(spikes, :unit);
            lf = groupby(theta,  [:tetrode]);
            sp = (G|>collect)[2]
            Threads.@threads for sp in G
                println("sp.unit[1] = ", sp.unit[1])
                lfkey = (;tetrode=cell_to_tet[sp.unit[1]])
                if lfkey ∉ keys(lf)|>collect.|>NamedTuple
                    continue
                end
                l = lf[lfkey]
                if isempty(lf)
                    continue
                end
                # println("extrema of lfp time:" ,    extrema(l.time).|>round, 
                #         "\nextrema of spike time:", extrema(sp.time).|>round)
                # Move phase over per tetrode
                DIutils.filtreg.register(l, sp; on="time", 
                                         transfer=[String(cycle)])
                cycles = convert(Vector{Union{Float32,Missing}}, 
                        sp[!,String(cycle)])
                ncycles = length(unique(cycles))
                if ncycles < 2
                    println("Only one-two cycle found for $(sp.unit[1])")
                else
                    println("Found $ncycles cycles for $(sp.unit[1])")
                end
                G[(;unit=sp.unit[1])][!,String(cycle)] .= cycles
            end
        end
        if added_lfp_field
            theta = theta[!, Not(String(cycle))]
        end
        GC.gc()
        if all(ismissing.(spikes[!, String(cycle)]))
            throw(ErrorException("No cycles found"))
        else
            spikes
        end
    end

    function isolated(spikes::AbstractDataFrame; kws...)
        prog = Progress(length(unique(spikes.unit)); 
                        desc="Adding isolation stats")
        func = x->(i=isolated(x; kws...);next!(prog);i)
        # combine(groupby(spikes, :unit), func)
        spikes = groupby(spikes, :unit)
        for cell in spikes
            try
                func(cell)
            catch exception
                print(exception)
                sleep(0.1)
            end
        end
        combine(spikes, identity)
    end
    isolated! = isolated

    """
        isolated(spikes::SubDataFrame; N=3, thresh=8, cycle_prop=:cycle, 
                 include_samples::Bool=false, overwrite=false)
    
    find isolated spikes in the manneer of Jai/Frank 2021
    
    # Arguments
    - `spikes::SubDataFrame`: dataframe of spikes
    - `N=3`: number of bins to look around
    - `thresh=8`: threshold for isolation
    - `cycle_prop=:cycle`: property to use for cycle
    - 1`include_samples::Bool=false`: whether to include the samples
    - `overwrite=false`: whether to overwrite existing columns
    
    # Returns
    - `SubDataFrame`: dataframe with isolation stats
    
    """
    function isolated(spikes::SubDataFrame; N::Int=3, thresh::Int=8,
    cycle::Symbol=:cycle, include_samples::Bool=false, overwrite=false)
        explore = setdiff(-N:N,0)
        println("Setting up unit $(spikes.unit[1])")
        if overwrite || !hasproperty(spikes, :isolated)
            spikes[!,:isolated] = Vector{Union{Missing,Bool}}(missing, size(spikes,1))
        end
        if overwrite || !hasproperty(spikes, :nearestcyc)
            spikes[!,:nearestcyc] = Vector{Union{Missing,Int32}}(missing, size(spikes,1))
        end
        if overwrite || !hasproperty(spikes, :meancyc)
            spikes[!,:meancyc] = Vector{Union{Missing,Float32}}(missing, size(spikes,1))
        end
        if overwrite || (include_samples && !hasproperty(spikes, :isosamples))
            spikes[!,:meancyc] = Vector{Union{Missing,Vector}}(missing, size(spikes,1))
        end
        spcycles = groupby(spikes, cycle)
        if all(ismissing.(spikes[!,cycle]))
            @warn "All cycles are missing" unit=spikes.unit[1]
        end
        #if length(cycles) > 1
            #@warn  "You only have 1 cycle"
        #end
        #(c,cycle) =  first(enumerate(cycles))
        print("Running unit ", spikes.unit[1])
        for (c,spcycle) in enumerate(spcycles)
            # find the N closest cycles
            explore_cycles = unique(max.(min.(c .+ explore, 
                                        [size(spcycles,1)]),[1]))
            if length(explore_cycles) < length(explore)
                spcycle.isolated .= false

                center_times = [abs(mean(c.time)-mean(spcycle.time)) 
                                for c in spcycles[explore_cycles]]
                order   = sortperm(center_times)
                nearest = explore_cycles[order] # TODO: see below
                cycle_prox = [abs(spcycle[1,cycle] - 
                              other_cyc[1,cycle])
                              for other_cyc in spcycles[nearest]]
                nearestcyc = @async minimum(cycle_prox)
                meancyc    = @async mean(cycle_prox)
                spcycle.meancyc    .= fetch(meancyc)
                cycle.nearestcyc .= fetch(nearestcyc)

            else
                center_times = [abs(mean(c.time)-mean(spcycle.time)) 
                for c in spcycles[explore_cycles]]
                order = sortperm(center_times)
                # TODO: make this match TODO above
                nearest = explore_cycles[order[1:N]] 
                cycle_prox = [abs(spcycle[1,cycle] - 
                                other_cyc[1,cycle])
                              for other_cyc in spcycles[nearest]]
                nearestcyc = @async minimum(cycle_prox)
                meancyc    = @async mean(cycle_prox)
                spcycle.meancyc    .= fetch(meancyc)
                spcycle.nearestcyc .= fetch(nearestcyc)
                isolated   = meancyc > thresh
                spcycle.isolated   .= isolated
                if include_samples
                    spcycle.isosamples .= [cycle_prox]
                end
            end
        end
        spikes
    end

    """
        rate_todataframe

    converts a rate dimarray (or count dimarray) to a dataframe, and populates columns labeling the
    times from other dataframes

    # Arguments
    - `X`: a DimArray with time as the first dimension
    - `registrant`: a tuple of (df, time_col, cols_to_transfer)
    """
    function rate_todataframe(X::DimArray, registrant::Tuple{DataFrame, String, Vector{String}}...)::DataFrame
        Table.from_dimarrary(X, registrant...)
    end

    """
        ripple_spikestats
    
    Add details about ripple events to spike data
    """
    function event_spikestats!(spikes::AbstractDataFrame,
                               events::AbstractDataFrame;
                eventname = "ripple")

        eventname = string(eventname)
        spikes[!, eventname] = Vector{Int32}(undef, size(spikes,1))
        spikes[!, eventname*"_phase"] = Vector{Float32}(undef, size(spikes,1))
        spikes[!, eventname*"_time"]  = Vector{Float32}(undef, size(spikes,1))

        event_of = DIutils.in_range_index(spikes.time, events; 
                                   start=:start, stop=:stop)
        event_of = convert(Vector{Int32}, event_of)
        spikes[!,eventname] = event_of

        spikes = subset(spikes, :ripple=>r->r .!= 0, view=true)

        # Add ripple phase
        spikes[:, eventname*"_start"] = 
            events[spikes[!,eventname], :start]
        spikes[:, eventname*"_stop"] = 
        events[spikes[!,eventname], :stop]
        spikes[!, eventname*"_phase"] = 
            (spikes.time .- spikes[!,eventname*"_start") ./ 
            (spikes[!,eventname*"_stop" .- spikes[!,eventname*"_start"])
        spikes[!, eventname*"_time"] = (spikes.time .- spikes[!,eventname*"_start"])

        return spikes
    end

    function event_spikestats!(spikes::GroupedDataFrame,
                               events::GroupedDataFrame;
             eventname = "ripple")
        spikes[!, eventname]          = Vector{Int32}(undef, size(spikes,1))
        spikes[!, eventname*"_phase"] = Vector{Float32}(undef, size(spikes,1))
        spikes[!, eventname*"_time"]  = Vector{Float32}(undef, size(spikes,1))
        K1 = keys(spikes) .|>  NamedTuple
        K2 = keys(events)  .|> NamedTuple
        K = intersect(K1,K2)
        for k in K
            spikes[k] = event_spikestats!(spikes[k], events[k], eventname=eventname)
        end
    end

end



