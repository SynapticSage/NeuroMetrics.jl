module task

    using DataFramesMeta
    """
        get_boundary(task; closed::Bool=false)
    acquires the boundary of a task
    # Input
    task    <- task dataframe
    closed  <- whether to close the boundary    
`   # Returns
    A boundary array, either closed or not
    """
    function get_boundary(task; closed::Bool=false)
        boundary = @subset(task, :name .== "boundary", :type .== "run")
        boundary = DataFrame(groupby(boundary, :epoch)[1])
        B = [boundary.x boundary.y]
        closed ? vcat(B, B[1,:]') : B
    end

end
