module timeshift

    __revise_mode__ = :evalassign
    using DataFrames
    import ..field
    export get_field_shift
    export field
    using ThreadSafeDicts
    using DataStructures
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    using Revise
    include("../table.jl")
    include("../shuffle.jl")
    using Infiltrator
    using Distributed
    using Dagger
    using Serialization
    using DrWatson
    using StatsPlots: @df
    using Plots
    using NaNStatistics
    using LoopVectorization

    export get_field_shift_shuffles, get_field_shift
    export to_dataframe, info_to_dataframe
    export fetch_best_fields
    export fetch_best_fields
    export func, correct, significant
    export crossval
    export ts_plotdir
    export plot_shifts
    export info_dataframe_and_cell_dataframe
    export save_mains, save_shuffles

    # -------------------- SHIFTING TYPES ---------------------------
    shift_func(data::DataFrame, shift::Real) = 
             transform(data, :time => (t->t.+shift) =>:time, copycols=false)
    const σ = shift_func

    # -------------------- SHIFTED Receptive Fields --------------------------
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; kws...)
        print("Using single")
        fieldobj = field.get_fields(σ(beh, shift), data; kws...)
    end

    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            multithread::Bool=true,
            postfunc::Union{Function,Nothing}=nothing,
            as::Type=OrderedDict,
            multi::Union{Bool, Symbol}=true,
            safe_dict::AbstractDict=ThreadSafeDict(),
            kws...)

        kws = (;dokde=false, kws...)

        if multi isa Bool
            multi = multi ? :thread : :single
        end
        @info "Starting multi=$multi"
        msg = "$multi timeshift-shuffles"

        p = Progress(length(shifts), desc="field shift calculations")
        @infiltrate
        if multi == :thread
            Threads.@threads for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = field.get_fields(σ(beh, shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        elseif multi == :single
            @showprogress 0.1 "shifts=$shifts" for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        elseif multi == :distributed
            @error "Not implemented"
        end
        safe_dict = Dict(safe_dict...)
        out = as(key=>pop!(safe_dict, key) for key in sort([keys(safe_dict)...]))
        return out
    end

    include("./timeshift/checkpoint.jl")
    include("./timeshift/dataframe.jl")
    include("./timeshift/operation.jl")
    include("./timeshift/plot.jl")
    include("./timeshift/shuffle.jl")

end
