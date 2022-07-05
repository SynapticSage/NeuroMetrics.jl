module shuffle

    import Field
    import Shuf
    using ThreadSafeDicts
    using DataFrames
    using DataStructures
    using Infiltrator
    using ProgressMeter
    using LoopVectorization
    using ..Timeshift: σ, shift_func

    export get_field_shift_shuffles

    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
    # COMPUTE THE a shuffle SETTING and distribute to downstream function
    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    function get_field_shift_shuffles(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
                preset::Union{Symbol, NamedTuple, AbstractDict},
                nShuffle::Int=100, 
                compute::Symbol=:single,
                postfunc::Union{Function,Nothing}=nothing,
                safe_dict::AbstractDict=ThreadSafeDict(),
                exfiltrateAfter::Real=Inf,
                get_field_kws...)::AbstractDict

        @info "Applying preset=$preset"
        initial_data_partials = Shuf.applyStandardShuffle(preset)
        _apply_partials(beh, data, shifts, initial_data_partials;
                        nShuffle, compute, postfunc,
                        safe_dict, exfiltrateAfter, get_field_kws...)
    end

    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
    # Functions for distributing partial functionals that generate data
    # from the user settings
    # ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
    function _apply_partials(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen,Vector{T}} where T <: Real,
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

        out = _run_partial_functional(beh, data, shifts, shuffle_data_generator;
                                      compute, nShuffle, postfunc, safe_dict,
                                      exfiltrateAfter, get_field_kws...)
    end

    function _apply_partials(beh::DataFrame, data::DataFrame,
                shifts::Union{StepRangeLen,Vector{T}} where T <: Real,
                initial_data_partials::T where T <: Function;
                nShuffle::Int=100, 
                compute::Symbol=:single,
                postfunc::Union{Function,Nothing}=nothing,
                safe_dict::AbstractDict=ThreadSafeDict(),
                exfiltrateAfter::Real=Inf,
                get_field_kws...)::AbstractDict

        partial = initial_data_partials
        shuffle_data_generator() = partial(data)

        out = _run_partial_functional(beh, data, shifts; 
                                      shuffle_data_generator, compute, nShuffle,
                                      postfunc, safe_dict, exfiltrateAfter,
                                      get_field_kws...)
    end

    function _run_partial_functional(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            shuffle_data_generator::Function,
            compute::Symbol, 
            nShuffle::Int,
            fieldfunc::Union{Function,Symbol}=Field.get_fields,
            postfunc::Union{Function,Nothing}=nothing,
            safe_dict::AbstractDict=ThreadSafeDict(),
            exfiltrateAfter::Real=Inf,
            get_field_kws...
        )

        fieldfunc = fieldfunc isa Symbol ? eval(fieldfunc) : fieldfunc

        # Collect sets we will iterate
        shuffle_shift_sets = collect(enumerate(Iterators.product(1:nShuffle, shifts)))

        function core_shuffle(shuffle::Int, shift::Real, P::Progress, i::Int)
            
            # Shuffle and run with shifted behavior
            data   = shuffle_data_generator()
            result = Field.get_fields(σ(beh,shift), data; get_field_kws...)

            # Apply post-processing
            if postfunc !== nothing
                result = postfunc(result)
            end

            # Push result and maybe exfiltrate
            push!(safe_dict, (;shift,shuffle)=>result)
            if mod(i-skips, exfiltrateAfter) == 0
                @info "chechpoint->exfiltrated"
                @exfiltrate
            end
            next!(P)
        end

        @info "Starting compute=$compute"
        msg = "$compute timeshift-shuffles"
        skips = 0
        P = Progress(length(shuffle_shift_sets), desc=msg)
        if compute == :thread
            Threads.@threads for (i, (shuffle,shift)) in shuffle_shift_sets
                if (;shift,shuffle) ∈ keys(safe_dict)
                    skips+=1
                    next!(P)
                    continue
                end
                core_shuffle(shuffle, shift, P, i)
            end
        elseif compute == :single
            @inbounds for (i,(shuffle,shift)) in shuffle_shift_sets
                if (;shift,shuffle) ∈ keys(safe_dict)
                    skips+=1
                    next!(P)
                    continue
                end
                core_shuffle(shuffle, shift, P, i)
            end
        else
            throw(ArgumentError("Unrecognized argument compute=$compute"))
        end
        safe_dict = Dict(safe_dict...)
        out = OrderedDict(key=>pop!(safe_dict, key) 
                          for key in sort([keys(safe_dict)...]))
    end

end
