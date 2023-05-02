module timeshift

    using DimensionalData
    using ...Timeshift.types
    using ...Timeshift.shiftmetrics
    using ...Field.metrics
    using DataFrames

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
    """
        usefulmetrics(M::DimArray, cells::AbstractDataFrame;
            filt_bad::Bool=false, brain_area::Union{String,Nothing}=nothing)

    Annotates shfited fields with useful metrics, and filters out bad cells
    or a specific brain area.
    # Arguments
    - `M::DimArray`: A matrix of shifted field objects
    - `cells::AbstractDataFrame`: A cell table
    - `filt_bad::Bool=false`: Whether to filter out bad cells
    - `brain_area::Union{String,Nothing}=nothing`: Whether to filter out cells
        from a specific brain area
    # Returns
    - `M::DimArray`: A matrix of shifted field objects with useful metrics
                     added to each shifted field object's dictionary
    """
    function usefulmetrics(M::DimArray, cells::AbstractDataFrame; 
        filt_bad::Bool=false, brain_area::Union{String,Nothing}=nothing,
        popunit::Bool=false)
        push_dims!(M)
        push_celltable!(M, cells)
        popunit ? pop_metric!(M, :unit) : nothing
        push_metric!(M, metrics.bitsperspike)
        push_metric!(M, metrics.totalcount)
        push_metric!(M, metrics.maxrate)
        push_shiftmetric!(M, best_tau!; metric=:bitsperspike)
        push_metric!(M, metrics.bitsperspikeold)
        push_shiftmetric!(M, best_tau!; metric=:bitsperspikeold)
        if filt_bad
            skcount=all(M[:totalcount] .> filt_requirements[:spikecount], 
                        dims=2)
            bps =   any(M[:bitsperspike] .> filt_requirements[:bitsperspike], 
                        dims=2)
            M = M[vec(skcount .&& bps) ,:]
        end
        if brain_area !== nothing
            M = M[vec(all(M[:area] .== brain_area, dims=2)), :]
        end
        M
    end
        

end
