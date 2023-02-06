module timeshift


    import Shuffle
    using DIutils
    import Timeshift: ShiftedField, ShiftedFields, AbsDictOfShiftOfUnit,
                      DictOfShiftOfUnit
    using Plots
    using Infiltrator
    using StatsPlots: @df
    using Statistics, NaNStatistics
    using DataFrames, DataFramesMeta
    using StatsBase
    using RecipesBase

    dosave = true

    #            -------------------------------------------------------
    #            --.--                                  o               
    #              |  ,   .,---.,---.    ,---.,---.,---..,---.,---.,---.
    #             |  |   ||   ||---'    |    |---'|    ||   ||---'`---.
    #             `  `---||---'`---'    `    `---'`---'`|---'`---'`---'
    #                `---'|                             |              
    #            -------------------------------------------------------
    #            Type recipes : recipes for timeshift types

    @userplot ShiftedFieldPlot
    @recipe function f(S::ShiftedFieldPlot)

        if length(S.args) == 1
            SF, shifts = S.args[1], nothing
        elseif length(S.args) == 2
            SF, shifts = S.args[1], S.args[2]
        else
            @error "Wrong argument number"
        end
        @assert SF isa ShiftedField
        @assert shifts isa Union{Vector, Nothing}

        shifts = shifts === nothing ? SF.keys : shifts
        D = SF.values

        legend := false
        layout := @layout [ grid(1, length(shifts)) ]
        grid := false
        aspect_ratio --> 1
        @info length(shifts)

        for (i,(s,d)) in enumerate(zip(shifts, D))
            @series begin
                subplot    := i
                seriestype := :heatmap
                c := :linear_kryw_5_100_c67_n256
                [d.grid.centers[1]...], [d.grid.centers[2]...], d.rate'
            end
        end
    end

    @userplot MetricShiftedFieldPlot
    @recipe function g(S::MetricShiftedFieldPlot; metric=nothing)

        if length(S.args) == 1
            SF, shifts = S.args[1], nothing
        elseif length(S.args) == 2
            SF, shifts = S.args[1], S.args[2]
        else
            @error "Wrong argument number"
        end
        @assert SF isa ShiftedField
        @assert shifts isa Union{Vector, Nothing}

        if metric === nothing
            #metric = unique()
        end

        D = SF.values
        M = SF.metrics
        shifts = shifts === nothing ? shifts.keys : shifts

        legend := false
        layout := @layout [ grid(1, length(shifts)) ]
        grid := false
        aspect_ratio --> 1
        @info length(shifts)

        for (i,(_,d)) in enumerate(zip(shifts, D))
            @series begin
                subplot    := i
                seriestype := :heatmap
                c := :linear_kryw_5_100_c67_n256
                [d.grid.centers[1]...], [d.grid.centers[2]...], d.rate'
            end
        end
    end

    # - 

end
