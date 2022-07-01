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

    # -----------
    # NEW SYSTEM|
    # -----------

    struct FixedGrid
        centers::Tuple
        edges::Tuple
        grid::Array
    end

    struct Occupancy
        count::Array{Float32}
        prob::Array{Float32}
    end

    struct Field
        count::Array{Float32}
        rate::Array{Float32}
        occ::Occupancy
        grid::FixedGrid
    end

    ReceptiveFields = AbstractDict{Field}

    # -----------
    # OLD SYSTEM|
    # -----------

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
