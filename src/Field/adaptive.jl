module adaptive

    using ..Field
    import ..Field: ReceptiveField
    import ..Field: get_boundary, resolution_to_width, return_vals
    import ..Field.metrics: MetricSet, push_metric!, pop_metric!
    import ..Field: metrics
    using DIutils, DIutils.binning
    import Filt

    using DataStructures
    using DataFrames
    import Load.utils: filterAndRegister, register
    import Base
    using LoopVectorization
    using ProgressLogging, ProgressMeter
    using Entropies: Probabilities
    using ProgressLogging
    using RecipesBase
    using Statistics
    using Polyester
    using ThreadSafeDicts
    using Infiltrator

    export GridAdaptive
    export AdaptiveOcc

    metric_def = [metrics.bitsperspike, metrics.totalcount, metrics.maxrate,
        metrics.maxcount, metrics.meanrate, metrics.coherence,
        metrics.argmax]

    using Plots
    using DataFrames: ColumnIndex
    CItype = Union{ColumnIndex,Vector{<:ColumnIndex}}
    CItype_plusNull = Union{ColumnIndex,Vector{<:ColumnIndex},Nothing}

    export yartsev

    #################################################
    ##### yartsev paper based  ####################
    #################################################
    thread_field_default = true
    thread_fields_default = false

    struct AdaptiveRF <: ReceptiveField
        grid::GridAdaptive
        occ::AdaptiveOcc
        count::Array{Int32}
        rate::Array{Float32}
        metrics::MetricSet
    end

    AdaptiveFieldDict = OrderedDict{<:NamedTuple,AdaptiveRF}
    #AdaptiveFieldDict(x) = OrderedDict{NamedTuple,   AdaptiveRF}(x)
    #AdaptiveFieldDict()  = OrderedDict{NamedTuple,   AdaptiveRF}()

    """
        yartsev(spikes, behavior, props; kws...)

    computes an adaptive grid and ratemap based on methods in yartsev papers
    """
    function yartsev(behavior::DataFrame, spikes::DataFrame, props::Vector;
        splitby::CItype_plusNull=[:unit],
        filters::Union{<:AbstractDict,Nothing}=nothing,
        metrics::Union{Function,Vector{Function},Nothing}=metric_def,
        thread_field::Bool=thread_field_default,
        thread_fields::Bool=thread_fields_default,
        prog_fields::Bool=false,
        grid_kws...)::Union{AdaptiveFieldDict,AdaptiveRF}
        #@info "yartsev" prog_fields
        if filters !== nothing
            if Filt.filters_use_precache(filters) &&
               Filt.missing_precache_output_cols(spikes, filters)
                @info "yartsev preacaching"
                spikes = Filt.precache(spikes, behavior, filters)
            end
            @info "filter and reg"
            @time behavior, spikes = filterAndRegister(behavior, spikes; filters,
                on="time", transfer=props,
                filter_skipmissingcols=true)
        else
            exists = intersect(Symbol.(props),
                propertynames(spikes))
            if length(exists) != length(props)
                behavior, spikes = register(behavior, spikes;
                    on="time", transfer=String.(props))

            end
        end
        @info "grid"
        @time grid = get_grid(behavior, props; grid_kws...)
        @info "occupancy"
        @time occ = get_occupancy(behavior, grid)
        @info "dropmissing"
        @time spikes = dropmissing(spikes, props)
        @info "fields"
        @time yartsev(spikes, grid, occ; splitby, metrics,
            prog_fields, thread_field, thread_fields,
            grid_kws...)
    end

    function yartsev(spikes::DataFrame, grid::GridAdaptive, occ::AdaptiveOcc;
        splitby::CItype=[:unit],
        metrics::Union{Function,Vector{Function},Nothing}=metric_def,
        thread_field::Bool=thread_field_default,
        thread_fields::Bool=thread_fields_default,
        prog_fields::Bool=false,
        grid_kws...)::Union{AdaptiveFieldDict,AdaptiveRF}

        #@info "yartsev" prog_fields
        get_adaptivefields(groupby(spikes, splitby), grid, occ;
            prog_fields, metrics, thread_field,
            thread_fields)
    end

    """
        get_adaptivefields(spikeGroups::GroupedDataFrame, props::Vector,
        grid::GridAdaptive; kws...)::AdaptiveFieldDict

    computes adaptive ratema based on a fixed grid derived from behavior
    """
    function get_adaptivefields(spikeGroups::GroupedDataFrame,
        grid::GridAdaptive, occ::AdaptiveOcc;
        thread_fields::Bool=thread_fields_default,
        prog_fields::Bool=false,
        metrics::Union{Function,Vector{Function},Nothing}=metric_def,
        kws...)::AdaptiveFieldDict
        keys_and_groups = collect(zip(Table.group.nt_keys(spikeGroups),
            spikeGroups))
        #@info "get_adaptivefields" prog_fields
        if thread_fields
            D = ThreadSafeDict{NamedTuple,AdaptiveRF}()
            Threads.@threads for (nt, group) in keys_and_groups
                D[nt] = get_adaptivefield(DataFrame(group), grid, occ; metrics, kws...)
            end
            D = OrderedDict(D)
        else
            D = OrderedDict{NamedTuple,AdaptiveRF}()
            prog = prog_fields ? Progress(length(keys_and_groups),
                desc="adaptive fields") : nothing
            for (nt, group) in keys_and_groups
                D[nt] = get_adaptivefield(DataFrame(group), grid, occ; metrics, kws...)
                prog_fields ? next!(prog) : nothing
            end
        end
        return D
    end

    """
        get_adaptivefield(X::DataFrame, props::Vector,
                          grid::GridAdaptive, occ::AdaptiveOc)::AdaptiveRF

    computes adaptive ratemap based on a fixed grid derived from behavior
    """
    function get_adaptivefield(spikes::DataFrame,
        grid::GridAdaptive, occ::AdaptiveOcc;
        metrics::Union{Function,Vector{Function},Nothing}=metric_def,
        thread_field::Bool=thread_field_default
    )::AdaptiveRF
        vals = Field.return_vals(spikes, grid.props)
        count = zeros(Int32, size(grid))
        if size(vals, 1) > 1
            if thread_field
                Threads.@threads for (index, (center, radius)) in collect(enumerate(grid))
                    if all((!).(isnan.(radius)))
                        @inbounds count[index] = sum(inside(vals, center, radius)) 
                    end
                end
            else
                @inbounds for (index, (center, radius)) in collect(enumerate(grid))
                    count[index] = sum(inside(vals, center, radius))
                end
            end
        end
        #@infiltrate all(count .== 0)
        count = reshape(count, size(grid))
        rate = @fastmath occ.camerarate * Float32.(count ./ occ.count)
        field = AdaptiveRF(grid, occ, count, rate, MetricSet())
        if metrics !== nothing
            for metric in metrics
                push_metric!(field, metric)
            end
        end
        field
    end



    ## --------
    ## UTILITIES
    ## --------
    function DIutils.dict.to_dict(F::AdaptiveRF)
        FF = Dict{Symbol,Any}()
        FF[:count] = F.count
        FF[:rate] = F.rate
        FF[:occ_prob] = reshape(F.occ.prob, size(F.grid.grid))
        FF[:occ_count] = F.occ.count
        FF[:grid_centers] = F.grid.centers
        FF[:props] = F.grid.props
        FF[:radii] = F.grid.radii
        FF
    end


    ## ------
    ## Skaggs
    ## ------
    #function skaggs_get_grid(behavior, props; width::Int, kws...)
    #    width = OrderedDict(prop=>width for prop in props)
    #    skaggs_get_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_get_grid(behavior, props; widths::Vector{<:Int}, kws...)
    #    width = OrderedDict(prop=>width for (prop,width) in zip(props,widths))
    #    skaggs_get_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_get_grid(behavior, props; width::OrderedDict, kws...)
    #    boundary = get_boundary(behavior, props)
    #    skaggs_get_grid(behavior, props; width, boundary, kws...)
    #end
    #function skaggs_get_grid(spikes, behavior, props; 
    #        thresh::Real=1, # Threshold in seconds
    #        dt=1/30, # Total time of sample
    #        radiusinc=0.1, # Spatial unit of RF
    #        width::OrderedDict, boundary::OrderedDict)
    #    vals_behavior = return_vals(behavior, props)
    #    vals_spikes   = return_vals(spikes, props)
    #    G = GridAdaptive(width, boundary, width)
    #    @avx for (index, center, radius) in cenumerate(G)
    #        while (sum(vals .< (center .+ radius)) * dt) < thresh
    #            radius += radiusinc
    #        end
    #        G.radii[index] = radius
    #    end
    #    G
    #end
    #"""
    #"""
    #function skaggs()
    #end
    
end
