module shuffle

    import ..Timeshift
    using ..Timeshift: σ, shift_func
    import Field
    import Field.preset: field_presets, return_preset_funcs
    import Field: adaptive
    import Shuf

    using ThreadSafeDicts
    using DataFrames
    using DataStructures
    using ProgressMeter
    using LoopVectorization
    using Infiltrator

    export get_field_shift_shuffles

    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
    # COMPUTE THE a shuffle SETTING and distribute to downstream function
    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    """
        get_field_shift_shuffles(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen, Vector{T}} where T <: Real; 
                fieldpreset::Union{Symbol, NamedTuple, AbstractDict},
                shufflepreset::Union{Symbol, NamedTuple, AbstractDict},
                nShuffle::Int=100, 
                compute::Symbol=:single,
                postfunc::Union{Function,Nothing}=nothing,
                safe_dict::AbstractDict=ThreadSafeDict(),
                exfiltrateAfter::Real=Inf,
                get_field_kws...)::AbstractDict

    Overarching function that dispatches first to find a data generation
    process. and then it dispatches that process and keywords to a final
    function that takes the data generator for the shuffle (be it a
    distribution or permutation) to a function that repeats generation and
    measurement for the number of shuffles requested.

    """
    function shifted_field_shuffles(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen, Vector{T}} where T <: Real,
                props::Vector; 
                shufflepreset::Union{Symbol, NamedTuple, AbstractDict},
                nShuffle::Int=100, 
                compute::Symbol=:single,
                postfunc::Union{Function,Nothing}=nothing,
                safe_dict::AbstractDict=ThreadSafeDict(),
                exfiltrateAfter::Real=Inf,
                get_field_kws...)::AbstractDict

        @info "Applying shufflepreset=$shufflepreset"
        initial_data_partials = Shuf.applyStandardShuffle(shufflepreset)
        _apply_partials(beh, data, shifts, props, initial_data_partials;
                        nShuffle, compute, postfunc,
                        safe_dict, exfiltrateAfter, get_field_kws...)
    end

    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
    # Functions for distributing partial functionals that generate data
    # from the user settings
    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
    """
    Apply the partial functions derived from shuffle presets to fill in a
    shuffle generated by some statistical DISTRIBUTION
    """
    function _apply_partials(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen,Vector{T}} where T <: Real,
                props::Vector,
                initial_data_partials::Tuple{<:Function,<:Function};
                nShuffle::Int=100, 
                compute::Symbol=:single,
                postfunc::Union{Function,Nothing}=nothing,
                safe_dict::AbstractDict=ThreadSafeDict(),
                exfiltrateAfter::Real=Inf,
                get_field_kws...)::AbstractDict

        partial, dist = initial_data_partials
        distribution  = dist(data)
        shuffle_data_generator() = partial(data, distribution)

        _run_partial_functional(beh, data, shifts, shuffle_data_generator;
                                      compute, nShuffle, postfunc, safe_dict,
                                      exfiltrateAfter, get_field_kws...)
    end

    """
    Apply the partial functions derived from shuffle presets to fill in a
    shuffle generated by a PERMUTATION
    """
    function _apply_partials(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen,Vector{T}} where T <: Real,
                props::Vector,
                initial_data_partials::T where T <: Function;
                nShuffle::Int=100, 
                compute::Symbol=:single,
                postfunc::Union{Function,Nothing}=nothing,
                safe_dict::AbstractDict=ThreadSafeDict(),
                exfiltrateAfter::Real=Inf,
                field_kws...)::AbstractDict

        partial = initial_data_partials
        shuffle_data_generator() = partial(data)

        _run_partial_functional(beh, data, shifts, props; 
                                      shuffle_data_generator, compute, nShuffle,
                                      postfunc, safe_dict, exfiltrateAfter,
                                      field_kws...)
    end

    """
    takes a datagenerator and field/shift options to compute shuffled and
    shifted fields.
    """
    function _run_partial_functional(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real,
            props::Vector; 
            shuffle_data_generator::Function,
            compute::Symbol, 
            nShuffle::Union{StepRangeLen, Int},
            safe_dict::AbstractDict=ThreadSafeDict{NamedTuple,Any}(),
            exfiltrateAfter::Real=Inf,
            skipproc::Bool=false,
            thread_field::Bool=true,
            precomputegridocc::Union{Bool,Nothing}=nothing,
            shiftbeh::Bool=false,
            field_kws...)


        # Collect sets we will iterate
        if nShuffle isa Int
            shuffle_sets = collect(enumerate(1:nShuffle))
        else
            shuffle_sets = collect(enumerate(nShuffle))
        end
        nShuffle = nShuffle isa Int ? (1:nShuffle) : nShuffle

        if precomputegridocc === nothing
            precomputegridocc = shiftbeh ? false : true
        end
        if precomputegridocc
            grid = adaptive.get_grid(beh, props; field_kws...) # TODO these would be different for fixed
            occ  = adaptive.get_occupancy(beh, grid)
            field_kws = (;field_kws..., grid, occ)
        end

        @showprogress "Shuffle" for s in nShuffle
            if skipproc && (;shuffle=s) ∈ keys(safe_dict)
                continue
            end
            data   = shuffle_data_generator()
            safe_dict[(;shuffle=s)] = Timeshift.shifted_fields(beh, data, shifts,
                                                       props;
                                                       overwrite_precache=true,
                                                       thread_field,
                                                       field_kws...)
            if mod(s, exfiltrateAfter)==0
                @exfiltrate
            end
        end

        safe_dict = Dict(safe_dict...)
        OrderedDict(key=>pop!(safe_dict, key) 
                          for key in sort([keys(safe_dict)...]))
    end

end
