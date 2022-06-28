module task

    using DataFramesMeta

    function get_boundary(task; closed::Bool=false)
        boundary = @subset(task, :name .== "boundary", :type .== "run")
        boundary = DataFrame(groupby(boundary, :epoch)[1])
        B = [boundary.x boundary.y]
        closed ? vcat(B, B[1,:]') : B
    end

end
