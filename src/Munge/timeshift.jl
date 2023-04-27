module timeshift

    using DimensionalData
    using ...Timeshift.types
    using ...Timeshift.shiftmetrics

    filt_requirements = Dict(
        :spikecount => 50,
        :bitsperspike => 0.5
    )

    function getshift(M::DimArray, ::Symbol)
        @assert(:bestshift_bitsperspike ∈ keys(M[1]),
        "Missing :bestshift_bitsperspike. Run Metrics.timeshift.besttau!")
        shifts = collect(M.dims[2])
        B = [findfirst(row) for row in 
            eachrow(M[:,1][:bestshift_bitsperspike].data .== shifts')]
        [unit[b] for (unit,b) in zip(eachrow(M),B)]
    end

    
    getshift(arrayOfFields::DimArray, s::T) where T <: Int = 
         arrayOfFields[:, arrayOfFields.dims[2].==s];


    function  usefulmetrics(SF::ShiftedFields)
        usefulmetrics(matrixform(SF))
    end
    function usefulmetrics(M::DimArray; 
        filt_bad::Bool=false, brain_area::Union{String,Nothing}=nothing)

        push_dims!(f)
        push_celltable!(f, cells)
        pop_metric!(f, :unit)
        push_metric!(f, metrics.bitsperspike)
        push_metric!(f, metrics.totalcount)
        push_metric!(f, metrics.maxrate)
        push_shiftmetric!(f, best_tau!; metric=:bitsperspike)
        push_metric!(f, metrics.bitsperspikeold)
        push_shiftmetric!(f, best_tau!; metric=:bitsperspikeold)
        if filt_bad
            skcount=all(f[:totalcount] .> filt_requirements[:spikecount], 
                        dims=2)
            bps = any(f[:bitsperspike] .> filt_requirements[:bitsperspike], 
                        dims=2)
            f = f[vec(skcount .&& bps) ,:]
        end
        if brain_area !== nothing
            f = f[vec(all(f[:area] .== brain_area, dims=2)), :]
        end
        f
    end
        

end
