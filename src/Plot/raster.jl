module raster
    #=
    #+NAME: raster.jl
    #+PURPOSE: Organizes my direct plot codes for rasters
    #           
    #           Section 1
    #           Plotting
    #
    # Links:
    # [[http://juliagraphics.github.io/Colors.jl/stable/namedcolors/][julia possible colors]]
    #
    # TODO option to separate dio and behavior into top level sepearte subplots
    =#

    using ColorSchemes, Colors, ColorSchemes, Gadfly
    using DataFrames, CSV, Statistics
    #using DrWatson
    using ProgressMeter, Debugger
    using Measures
    import DIutils.Table as Table

    __revise_mode__ = :eval

    StringSymbol = Union{String, Symbol}


    norm(x, y) = (maximum(y)-minimum(y)).*(x.-minimum(x))./(maximum(x)-minimum(x)).-minimum(y);
    dio_colors = repeat([colorant"white"], 5)
    padding = 0.05;
    skipnan(x) = Iterators.filter(!isnan, x)
    nanmax(x) = maximum(skipnan(x))
    mm = Measures.mm
    cm = Measures.cm

    function _setup_gadfly(spikes, beh,
            timegroup, timegroup_count;
            sortY=nothing,
            sortColor=nothing,
            behavior_vars=nothing, 
            behavior_geom=Geom.line,
            norm_behvar::Bool=false,
            colorTimeGroup=nothing,
            dio_vars=nothing,
        )

        DEBUG = Dict()
        if timegroup_count isa String
            timegroup_count = parse(Int, timegroup_count);
        end
        if sortY == nothing
            sortY = :unit;
        end
        if sortColor == nothing
            sortColor = :unit;
        end

        if timegroup_count != -1
            spikes = Table.select_row_group(spikes, timegroup, timegroup_count);
            beh    = Table.select_row_group(beh,    timegroup, timegroup_count);
        end
        if isempty(spikes)
            throw(InvalidStateException("Empty spikes", :spikes))
        elseif isempty(beh)
            throw(InvalidStateException("Empty beh", :beh))
        end
        if colorTimeGroup != nothing
            trajs = Table.get_periods(spikes, colorTimeGroup);
            runTrajLayer = layer(trajs, xmin=:start, 
                                 xmax=:end, color=:traj, alpha=[0.35], 
                                 Geom.vband);
            trajLayerElements = [Scale.color_discrete(), runTrajLayer];
        else
            trajLayerElements = [];
        end
        nonmissing = (!).(ismissing.(beh[!,behavior_vars]));
        if norm_behvar
            beh[nonmissing,behavior_vars] = norm(beh[nonmissing, behavior_vars], spikes[!,sortY])
        end
        behaviorLayer = Gadfly.layer(beh[nonmissing,:], 
                                    x=:time, y=behavior_vars, 
                                color=[colorant"red"], alpha=[0.5],
                                Gadfly.style(line_width=1*Measures.mm),
                                behavior_geom);

        # Plot DIO
        if dio_vars ≠ nothing
            maxY = maximum(spikes[!,sortY]);
            if dio_vars isa String
                if dio_vars in names(beh)
                    dio_vars = [dio_vars];
                else
                    dio_vars = Table.vec_of_matching_colnames(beh, dio_vars);
                end
            end
            dioLayer = Vector{Vector{Gadfly.Layer}}()
            DEBUG["dio_vars"] = dio_vars
            DEBUG["maxY"] = maxY
            for i = 1:length(dio_vars)
                dio   = replace(beh[!, dio_vars[i]] * maxY * (1+padding), missing=>NaN, 0=>NaN);
                alpha = replace(beh[!, dio_vars[i]], missing=>0);
                if any(beh[!,dio_vars[i]] .== 1)
                    println("diovar $i $(nanmax(beh[!,dio_vars[i]]))")
                    DEBUG["dio_$i"]= ((beh[!,:time], dio))
                    push!(dioLayer, layer(beh, x=:time, y=dio, alpha=alpha, color=[dio_colors[i]], Geom.line))
                end
            end
        else
            dioLayer = [];
        end

        spike_tick_style = Gadfly.style(point_shapes=[Gadfly.Shape.vline], 
                            default_color=colorant"white",
                            line_width=0mm,
                            point_size=0.5mm);
        return (spikes, trajLayerElements, behaviorLayer, spike_tick_style, dioLayer, DEBUG)
    end

    function _set_rank(rank=nothing, rankY=nothing, rankColor=nothing)
        if rank == "none"
            rank = false
        end
        if rankY == "none"
            rankY = false
        end
        if rankColor == "none"
            rankColor = false
        end
        if rankColor == nothing
            rankColor = false
        end
        if rankY == nothing
            rankY = rank;
        end
        if rankColor == nothing
            rankColor = rank;
        end
        (rankY, rankColor)
    end

    function _unit_specific_column(spikes, column, rank)
        if !(column in names(spikes))
            colnames = Table.vec_of_matching_colnames(spikes, column)
            if isempty(colnames)
                throw(ArgumentError("No matching to match=$column"))
            end
            if rank
                colnames  = [colname for colname in colnames if
                             occursin("_argcenter",colname)]
            else
                colnames  = [colname for colname in colnames if
                             occursin("_center",colname)]
            end
        else
            if rank
                colnames  = column * "_argcenter"
            else
                colnames  = column * "_center"
            end
        end
        return colnames
    end

    function gadfly_subplot(spikes::DataFrame, beh::DataFrame,
            timegroup::String, timegroup_count::Int;
            sortY::Union{Nothing,String}=nothing,
            sortColor::Union{Nothing,String}=nothing,
            subplot::Union{Nothing,String}=nothing,
            group_tabledim::Symbol=:row,
            stackdir::String="v",
            rankOrder::Union{Bool,Nothing,String}=nothing,
            rankY::Union{Bool,Nothing,String}=nothing,
            rankColor::Union{Bool,Nothing,String}=nothing,
            debug::Bool=false,
            norm_behvar::Bool=false,
            behavior_geom=Geom.line,
            dio_vars="poke",
            other...)

        DEBUG = Dict()

                                             
        #                       |              
        #             ,---.,---.|--- .   .,---.
        #             `---.|---'|    |   ||   |
        #             `---'`---'`---'`---'|---'
        #                                 |    

        if sortY == nothing; sortY = :unit; end
        if sortColor == nothing; sortColor = :unit; end
        rY, rColor = raster._set_rank(rankOrder, rankY, rankColor);
        sY, sColor = raster._unit_specific_column(spikes, sortY, rY),
                     raster._unit_specific_column(spikes, sortColor, rColor)

        # TIME GROUPING
        if timegroup_count isa String
            timegroup_count = parse(Int, timegroup_count);
        end
        if timegroup_count != -1
            S = Table.select_row_group(spikes, timegroup, timegroup_count);
            B    = Table.select_row_group(beh,    timegroup, timegroup_count);
        end
        if isempty(spikes)
            throw(InvalidStateException("Empty spikes", :spikes))
        elseif isempty(beh)
            throw(InvalidStateException("Empty beh", :beh))
        end

        # Color time group?
        if colorTimeGroup != nothing
            trajs = Table.get_periods(spikes, colorTimeGroup);
            runTrajLayer = layer(trajs, xmin=:start, 
                                 xmax=:end, color=:traj, alpha=[0.35], 
                                 Geom.vband);
            trajLayerElements = [Scale.color_discrete(), runTrajLayer];
        else
            trajLayerElements = [];
        end

        # Behavior layer
        nonmissing = (!).(ismissing.(B[!,behavior_vars]));
        if norm_behvar
            B[nonmissing,behavior_vars] = norm(B[nonmissing, behavior_vars], spikes[!,sortY])
        end
        behaviorLayer = Gadfly.layer(B[nonmissing,:], 
                                    x=:time, y=behavior_vars, 
                                color=[colorant"red"], alpha=[0.5],
                                Gadfly.style(line_width=1mm),
                                behavior_geom);

        # Plot DIO
        if dio_vars ≠ nothing
            maxY = maximum(spikes[!,sortY]);
            if dio_vars isa String
                if dio_vars in names(B)
                    dio_vars = [dio_vars];
                else
                    dio_vars = Table.vec_of_matching_colnames(B, dio_vars);
                end
            end
            dioLayers = Vector{Vector{Gadfly.Layer}}()
            DEBUG["dio_vars"] = dio_vars
            DEBUG["maxY"] = maxY
            for i = 1:length(dio_vars)
                dio   = replace(B[!, dio_vars[i]] * maxY * (1+padding),
                                missing=>NaN, 0=>NaN);
                alpha = replace(B[!, dio_vars[i]], missing=>0);
                if any(B[!,dio_vars[i]] .== 1)
                    print(i)
                    push!(dioLayers, layer(B, x=:time, y=dio, alpha=alpha, 
                                           Gadfly.style(line_width=4mm),
                                           color=[dio_colors[i]], Geom.line))
                    #b = B[1,:]
                    #push!(dioLayers, 
                    #      layer(x=[b.time], y=[dio[(!isnan).(dio)][1]], label=[dio_vars[i]],
                    #            size=[8pt], Geom.label)
                    #     )
                end
            end
        else
            dioLayers = [];
        end

        spike_tick_style = Gadfly.style(point_shapes=[Gadfly.Shape.vline], 
                            default_color=colorant"white",
                            line_width=0mm,
                            point_size=0.5mm);


        if subplot == nothing
            spikes_sep = [S];
        else
            spikes_sep = groupby(S, subplot);
        end

        multiysort = sY isa Vector && length(sY) > 1
        if multiysort
            println("More than 1")
            println(sY)
            iterate = Iterators.product(1:length(spikes_sep), enumerate(sY))
            ylen = length(sY)
        else
            print("Not more than 1")
            iterate = 1:length(spikes_sep);
            ylen = 1
        end
        println(collect(iterate))

        plots = repeat([plot()], length(spikes_sep), ylen);
        for unpack in iterate

            if multiysort
                i, (y,ysort) = unpack
            else
                i = unpack
                ysort = sY
                y=1
            end

            dat = DataFrame(spikes_sep[i])
            dat = dat[(!isnan).(dat[!,ysort]),:]
            if any(size(dat) .== 0)
                print("Zero")
            end

            neuralLayer = layer(dat, color=sColor, alpha=[0.5],
                                Gadfly.style(highlight_width=0mm),
                                spike_tick_style, Geom.point, y=ysort, x=:time); 

            thing = Gadfly.plot( neuralLayer, behaviorLayer, dioLayers...)
            plots[i,y] = thing;
        end

        # Stacking
        if stackdir == "v"
            plots = vstack(plots...)
        elseif stackdir == "h"
            plots = hstack(plots...)
        else
            throw(ArgumentError("stackdir cannot be $stackdir"));
        end

        if debug
            return plots, DEBUG
        else
            return plots
        end

    end

    # from [https://github.com/GiovineItalia/Gadfly.jl/issues/1169]
    function vstack2(plots::Vector{Plot}; spacing::Float64=0.0, heights::Vector{<:Real}=Float64[])
        n = length(plots)
        heights==Float64[] && (heights = fill(1/n, n))
        sumh = cumsum([0;heights])
        vpos = sumh + spacing.*[0; 1:n-1; n-1]
        M = [(context(0,v,1,h), render(p)) for (v,h,p) in zip(vpos[1:n], heights, plots)]
        return compose(context(units=UnitBox(0,0,1,vpos[n+1])), M...)
    end

    function set_default_theme(;dark=false)
        if dark
            extra = (Gadfly.dark_theme,)
        else
            extra = ()
        end
        bigger_fonts = Theme(extra...,
              highlight_width=0pt,
              key_position=:none,
              major_label_font_size=18pt,
              minor_label_font_size=14pt,
              key_title_font_size=18pt,
              key_label_font_size=14pt);
        Gadfly.push_theme(bigger_fonts);
        nothing
    end

    function populate_sort_fields(spikes::DataFrame, beh::DataFrame,
            props::Vector{String})
        spikes = dropmissing(spikes);
        @assert size(spikes,1) != 0
        # add sortable properties for cells
        println("Adding cell sort properties")
        spikes = Table.add_sort_properties(spikes, beh, props);
        areawise = groupby(spikes, "area");
        A = Vector{DataFrame}([])
        # add sortable properties for cells by area
        println("Adding cell sort properties by area")
        @showprogress "Iterating areas" for a = 1:length(areawise)
            A = [A; Table.add_sort_properties(DataFrame(areawise[a]), beh, props; 
                                              modifier="area_")];
        end
        spikes = vcat(A...)
        return spikes
    end

    function preproc_toggle_subplotSpecOrder(subplotToggle, metricToggle, sortby, colorby)
        # if splitby, switch our ordering variable, UNLESS, not arg
        if subplotToggle in ["sortby", "both"] && xgroup != nothing && !(metricToggle in ["sortby","both"])
            sortby = xgroup * "_" * sortby;
        elseif subplotToggle in ["colorby", "both"] && xgroup != nothing && !(metricToggle in ["colorby","both"])
            colorby = xgroup * "_" * colorby;
        end
        return sortby, colorby
    end

    function preproc_toggle_metricRank(toggle, sortby, colorby)
        # if splitby, switch our ordering variable, UNLESS, not arg
        if metricToggle in ["sortby", "both"] && xgroup != nothing
            replace!(sortby,"argcenter"=>"center")
        elseif metricToggle in ["colorby", "both"] && xgroup != nothing
            replace!(colorby,"argcenter"=>"center")
        end
        return sortby, colorby
    end

    #function gadfly_subplot(raster, beh,
    #        timegroup, timegroup_count;
    #        sortby=nothing,
    #        colorby=nothing,
    #        xgroup=nothing,
    #        ygroup=nothing,
    #        other...)
    #    #=
    #    DOESN"T WORK
    #    =#
    #    (raster,
    #     trajLayerElements,
    #     behaviorLayer,
    #     spike_tick_style) =_setup_gadfly(raster, beh, timegroup, timegroup_count; other...)
    #    neuralLayer = layer(Geom.point, y=sortby, x=:time, );
    #    subplotGrid = Geom.subplot_grid(neuralLayer, behaviorLayer);
    #    if xgroup != nothing && ygroup != nothing
    #        newview = plot(raster, spike_tick_style, 
    #                       xgroup=xgroup, ygroup=ygroup,
    #                       y=sortby, x=:time, 
    #                       color=colorby, alpha=[0.5], 
    #                       subplotGrid
    #                  );
    #    elseif xgroup != nothing && ygroup == nothing
    #        newview = plot(raster, spike_tick_style, 
    #                       xgroup=xgroup,
    #                       y=sortby, x=:time, 
    #                       color=colorby, alpha=[0.5], 
    #                       subplotGrid,
    #                  );
    #    elseif xgroup == nothing && ygroup != nothing
    #        newview = plot(raster, spike_tick_style, 
    #                       ygroup=ygroup,
    #                       y=sortby, x=:time, 
    #                       color=colorby, alpha=[0.5], 
    #                       subplotGrid,
    #                  );
    #    elseif xgroup == nothing && ygroup == nothing
    #        newview = plot(raster, spike_tick_style, 
    #                       y=sortby, x=:time, 
    #                       color=colorby, alpha=[0.5], 
    #                       subplotGrid,
    #                  );
    #    end
    #end

    #function gadfly(raster, beh,
    #        timegroup, timegroup_count;
    #        sortby=nothing,
    #        colorby=nothing,
    #        other...
    #    )
    #    (raster,
    #     trajLayerElements,
    #     behaviorLayer,
    #     spike_tick_style) =_setup_gadfly(raster, beh, timegroup, timegroup_count; other...)
    #
    #    coordinate = Coord.cartesian(xmin=minimum(raster.time), 
    #                                 xmax=maximum(raster.time), 
    #                                 ymax=maximum(raster[:,sortby]));
    #    newview = plot(raster, behaviorLayer, coordinate,
    #         x=:time, y=sortby, color=colorby, alpha=[0.5], spike_tick_style,
    #         Geom.point);
    #    newview
    #end

    function plotjl()
    end

end
