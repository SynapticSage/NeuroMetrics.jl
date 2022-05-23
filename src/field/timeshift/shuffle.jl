
"""
shuffled version of get_field_shift
"""
function get_field_shift_shuffles(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            multi::Union{Symbol,Bool}=true,
            nShuffle::Int=100, 
            shuffle_pos::Union{Tuple,NamedTuple}=(;),
            shuffle_kws::NamedTuple=(;),
            shuffle_func=shuffle.by,
            postfunc::Union{Function,Nothing}=nothing,
            safe_dict::AbstractDict=ThreadSafeDict(),
            exfiltrateAfter::Real=Inf,
            kws...)::AbstractDict

    kws = (;dokde=false, kws...)

    sets = collect( Iterators.enumerate(Iterators.product(1:nShuffle,
                                                          shifts)))

    # Generate the distribution functoin for all shuffles (expensive to do for each)
    if isempty(shuffle_pos) || !(shuffle_pos[1] isa Distribution)
        shuffle_kws = (; shuffledist_df=beh, shuffle_kws...)
        distribution = shuffle._create_distribution(data, 
                                                    shuffle_pos...;
                                                    shuffle_kws...)
        shuffle_pos = (;shuffle_pos..., distribution)
        @info "distribution=$distribution"
    end
    
    if multi isa Bool
        multi = multi ? :thread : :single
    end

    @info "Starting multi=$multi"
    msg = "$multi timeshift-shuffles"
    skips = 0
    P = Progress(length(sets), desc=msg)
    if multi == :thread
        Threads.@threads @inbounds for (i, (shuffle,shift)) in sets
            if (;shift,shuffle) ∈ keys(safe_dict)
                skips+=1
                continue
            end
            data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
            result = field.get_fields(σ(beh,shift), data; kws...)
            if postfunc != nothing
                result = postfunc(result)
            end
            push!(safe_dict, (;shift,shuffle)=>result)
            next!(P)
            if mod(i-skips, exfiltrateAfter) == 0
                @info "chechpoint->exfiltrated"
                @exfiltrate
            end
        end
    elseif multi == :distributed
        @assert nprocs() > 1 msg
        @showprogress 0.1 msg for @inbounds (i,(shuffle,shift)) in sets
            data = Dagger.spawn(shuffle_func, data, shuffle_pos...; shuffle_kws...)
            @debug "Dagger 1"
            result = Dagger.spawn(field.get_fields, σ(beh,shift), data; kws...)
            @debug "Dagger 2"
            if postfunc != nothing
                result = Dagger.@spawn postfunc(result)
                @debug "Dagger 3"
            end
            push!(safe_dict, (;shift,shuffle)=>result)
        end
    elseif multi == :single
        @inbounds for (i,(shuffle,shift)) in sets
            if (;shift,shuffle) ∈ keys(safe_dict)
                skips+=1
                next!(P)
                continue
            else
                #@info (; i,shift,shuffle)
            end
            data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
            result = field.get_fields(σ(beh,shift), data; kws...)
            if postfunc != nothing
                result = postfunc(result)
            end
            push!(safe_dict, (;shift,shuffle)=>result)
            if mod(i-skips, exfiltrateAfter) == 0
                @info "chechpoint->exfiltrated"
                @exfiltrate
            end
            next!(P)
        end
    else
        throw(ArgumentError("Unrecognized argument multi=$multi"))
    end
    safe_dict = Dict(safe_dict...)
    out = OrderedDict(key=>pop!(safe_dict, key) 
                      for key in sort([keys(safe_dict)...]))
    return out

end
