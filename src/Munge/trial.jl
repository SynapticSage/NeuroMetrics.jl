"""
    trial

trial-based transformations


outputs https://github.com/JuliaArrays/AxisArrays.jl types
"""
module trial

    SymStr = Union{Symbol, String}

    """
        tensor_pointproc

    create a tensor of a point process arrayed with some set of N-dimensions
    """
    function tensor_pointproc(X::DataFrame)
    end

    """
        tensor_pointproc

    create a tensor of a continuous variable from dataframe form

    ### Input
    
    `X`         -- DataFrame we'll pull a tensor out of
    `dims`      -- Columns that will become the dimension of the tensor
    `var`       -- Column that will become the measurement
    `equalize`  -- (optional, default=false) List of columns that we will force to be equal in length
    `equalmeth` -- (optional, default=:first) Method of equalization

    ### Output
    `T` -- Tensor; either a square array or stacked arrays of unequal dim
     
    """
    function tensor_continuous(X::DataFrame, dims::Vector, var::<:SymStr;
            equalize::Bool=false, equalmeth::Symbol=:first)
        
    end

    function equalize
    end

end
