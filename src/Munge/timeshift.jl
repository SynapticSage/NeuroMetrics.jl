module timeshift

    using DimensionalData

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

end
