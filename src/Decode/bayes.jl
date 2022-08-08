"""
    bayes

Tools for implementing bayesian decoding using Fields.jl and 
Timeshift.jl elements.
"""
module bayes

    using StatsBase
    using Entropies
    using AxisArrays
    using ProgressMeter

    import Timeshift

    function posterior(spikes::DataFrame, fields::ShiftedFields; binsize=0.02)

        fields = Timeshift.types.matrixform(fields)
        occ = fields[1,1].occ
        prob = reshape(occ.prob, size(occ.count))
        fields = Timeshift.types.tensorform(fields)
        counts = Munge.spiking.tocount(spikes; binsize)

        posterior(counts, prob, block; binsize)

    end
    function posterior(counts::AxisArrray, prob::AbstractArray, fields::AbstractArray;
        binsize=0.001)

        axs = Dict(zip(axisnames(fields), fields.axes))
        countaxes = Dict(zip(axisnames(counts), counts.axes))
        axs[:unit] = Axis{:unit}(intersect(axs[:unit], countaxes[:unit]))

        if :shift in keys(axs)
            size_prob = collect(size(prob))
            size_prob .= 1
            prob = repeat(prob, size_prob..., length(axs[:shift]))
            prob = permutedims(prob, (3, 1, 2))
            prob ./= length(axs[:shift])
        end

        # TODO ordrering possibly not correct
        F = [fields[unit=(axs[:unit] .== n)][1,:,:,:] for n in axs[:unit]]
        C = [counts[unit=(axs[:unit] .== n)][:,1,] for n in axs[:unit]]


        return prob
    end

    function tofloat16(F::Vector{<:AxisArray})::Vector{AxisArray}
        Q = Vector{AxisArray}(undef, length(F))
        for i in 1:length(F)
            @info i
            Q[i] = AxisArray(convert(Array{Float16}, F[i].data), F[i].axes)
        end
        Q
    end

    function decode_vectorized_log(prob, F::Vector, C::Vector; binsize)
        prog=Progress(length(F); desc="Bayesian decode")
        prob = log.(prob)
        prob .+= (-binsize .* sum(F)) # multiply prior by exponential term
        tmp = reduce(sum,
                     begin
                         F=field_raised_to_count_chunked(f, c, .+; chunk=100000)
                         next!(prog)
                         F
                     end
                         for (f,c) in zip(F, C))
        prob .+= tmp   # then by product of field ter

        prob = exp.(log)
        prob ./= sum(prob, dims=(collect(2:ndims(prob))...)) # normalize distros to 1
    end

    function decode_vectorized(prob, F::Vector, C::Vector; binsize)
        prog=Progress(length(F); desc="Bayesian decode")
        prob .*= exp.(-binsize .* sum(F)) # multiply prior by exponential term
        tmp = reduce(prod,
                     begin
                         F=field_raised_to_count_chunked(f,c; chunk=100000)
                         next!(prog)
                         F
                     end for (f,c) in zip(F, C))
        prob .*= tmp   # then by product of field ter
        prob ./= sum(prob, dims=(collect(2:ndims(prob))...)) # normalize distros to 1
    end

    function field_raised_to_count(f::Union{Array,AxisArray}, 
                                   c::Union{Array,AxisArray}, op=.^)
        op(f, permutedims(c[:,:,:,:], (2,3,4,1)))
    end

    function field_raised_to_count_chunked(f::Union{Array,AxisArray}, 
                                   c::Union{Array,AxisArray}, op=.^;
                                   chunk=100000)
        P = []
        c = permutedims(c[:,:,:,:], (2,3,4,1))
        @showprogress for t in 1:chunk:length(c)
            append!(P, op(f, c[:,:,:,t:chunk-1]))
        end
    end

    function bin()
        #function h2d(thing::DataFrame, props::Vector{String}; grid=(),
        #        hist2dkws=Dict())
        #    thing = dropmissing(thing);
        #    P = tuple([convert(Vector{Float64}, x) for x in
        #               eachcol(thing[:, props])]...)
        #    return fit(Histogram, P, grid; hist2dkws...)
        #end
    end

end
