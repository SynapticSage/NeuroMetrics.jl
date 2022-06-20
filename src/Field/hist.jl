module hist

    using ..Field
    import Utils

    using DataFrames
    using StatsBase: fit
    using StatsBase
    using Infiltrator

    function h2d(thing::DataFrame, props::Vector{String}; grid=(),
            hist2dkws=Dict())
        thing = dropmissing(thing);
        P = tuple([convert(Vector{Float64}, x) for x in eachcol(thing[:, props])]...)
        return fit(Histogram, P, grid; hist2dkws...)
    end

    function fields(data::DataFrame, beh::DataFrame;
            splitby::Union{Nothing,Vector{Symbol},Vector{String},String}="unit",
            resolution::Union{Vector{Int},Int}=50, 
            savemem::Bool=true,
            props::Vector{String}=["x","y"],
            gaussian::Real=0.0
        )
        if gaussian isa Int
            gaussian = convert(Float64, gaussian);
        end

        grid = Field.getSettings(beh, props, 
                                resolution=resolution,
                                settingType="hist");
        @debug "Grid=$grid"
        @debug "resolution=$resolution"

        behDist = Field.hist.h2d(beh, props, grid=grid);
        
        ith_group(i) = DataFrame(groups[i])
        field_of_group(i) = hist.h2d(ith_group(i), props, grid=grid)
        if splitby isa Vector{Symbol}
            splitby = String.(splitby)
        end
        if all(in.(splitby, [names(data)]))
            groups = groupby(data, splitby)
            ugroups =
            collect(zip(sort(collect(groups.keymap),by=x->x[2])...))[1]
            ugroups = 
            [(;zip(groups.cols,ugroup)...) for ugroup in ugroups]
            #spikeDist = Dict(ugroups[i] => field_of_group(i) 
            #                 for i ∈ 1:length(groups))
            spikeDist = Dict{typeof(ugroups[1]),Any}();
            @inbounds for i ∈ 1:length(groups)
                if savemem
                    spikeDist[ugroups[i]] = Float32.(field_of_group(i).weights)
                else
                    spikeDist[ugroups[i]] = field_of_group(i).weights
                end
            end
        else
            if splitby != nothing
                @warn "splitby not empty, yet not all columns in df"
            end
            spikeDist = hist.h2d(data, props, grid=grid)
        end

        # Transform behavior
        behDist_nanWhere0 = copy(behDist.weights);
        if savemem
            behDist_nanWhere0 = convert.(Float32, behDist_nanWhere0);
        else
            behDist_nanWhere0 = convert.(Float64, behDist_nanWhere0);
        end
        behzeroinds = behDist_nanWhere0 .== 0
        behDist_nanWhere0[behzeroinds] .= NaN;
        behDist.weights[behDist.weights .== 0] .= 1;
        
        return (hist=spikeDist, grid=grid, occ=behDist_nanWhere0,
                occzeroinds=behzeroinds)
    end
end
