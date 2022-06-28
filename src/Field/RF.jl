"""
RF

Receptive Field module

Module for typing field calculations

Work-in-progress

1. method to convert current field calcs into RF type
    1. to_RF()
    2. from_RF()
2. modules that have to adapt to the type
    1. info
    2. model
    3. operation
    4. plot
    5. recon/recon_process
"""
module RF

    struct ReceptiveField
            Cₕ::Array{Float32}
            Cₖ::Union{Array{Float32},Nothing}
            occ::Array{Float32}
            occR::Array{Float32, Nothing}
            occzeroinds::Array{Float32}
            cgrid::Tuple
            egrid::Tuple
            gridh::Tuple
            gridk::Tuple
            dims::Union{Vector{String}, Vector{Symbol}}
    end

    const ReceptiveFields = AbstractDict{ReceptiveField}

    struct Old_ReceptiveFields
            Cₕ::Dict
            Cₖ::Union{Dict, Nothing}
            occ::Array{Float32}
            occR::Array{Float32, Nothing}
            occzeroinds::Array{Float32}
            cgrid::Tuple
            egrid::Tuple
            gridh::Tuple
            gridk::Tuple
            dims::Union{Vector{String}, Vector{Symbol}}
    end

    _dictfields =  [:Cₕ, :Cₖ, :Rₕ, :Rₖ]

end
