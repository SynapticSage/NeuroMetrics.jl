"""
`plot`

# proposed structure

`show_fieldgroups` : plot a group of field collections, who are indexed by group keys

`show_fields` : plot a collection of fields indexed by keys

`show_field`  : plot a single field and its key

"""
module plot

    using ..field
    import ..utils
    using Plots, LaTeXStrings, Measures
    using Statistics
    using ProgressMeter

    function transform_key(key; splitter::String="\n")
        tk(key) = replace(string(key), "("=>"",",)"=>"",", "=>splitter,
                                     ")"=>"", "\""=>"", " "=>"")
        if key isa Nothing
            annotation = ""
        elseif key isa String
            annotation = key
        else
            annotation = tk(key)
        end
    end

    function show_fieldgroups(group; groupgrid="row", as=Plots.plot,
            show_field_kws...)
        if groupgrid == nothing
            groupgrid = "rowcol"
        end
        if groupgrid isa String
            key = Tuple(keys(group))[1]
            N = length(group)
            if groupgrid == "row"
                groupgrid = (N, 1)
            elseif groupgrid == "col"
                groupgrid = (1, N)
            elseif groupgrid == "rowcol"
                groupgrid = (ceil(sqrt(N)), floor(sqrt(N)))
            end
        end
        keys_ = Tuple(keys(group))
        N = length(keys_)
        if as == Plots.plot
            extrakws = (;nplots_factor=N, show_field_kws...)
        else
            extrakws = (; show_field_kws...,)
        end
        plots = [show_fields(group[key]; keyappend=key, extrakws...)
                 for key in keys_]
        if as == Plots.plot
            grid = [p.layout for p in plots]
            @debug "show_fieldgroups: gridsize=$(size(grid))"
            grid = reshape(grid, groupgrid)
            plots=Plots.plot(plots..., grid=grid)
        elseif as == Dict || as == NamedTuple
            plots = as(zip(keys_,plots))
        else
            plots = as(plots)
        end
        return plots
    end
    function save_dict_fieldplots(plots, name; ext="pdf")
        @showprogress for (key, plot) in plots
            keyname = name * "_" * transform_key(key, splitter=",")
            if ext isa String
                ext = [ext]
            end
            for e in ext
                savefig(plot, keyname * ".$e")
            end
        end
    end
    function selectbackground()
        ndim = ndims(operation.selectind(F, 1))
    end
    function show_fields(F::Dict; fontscale=true, background=:grey30,
            textcolor=:white, as::Union{Type,<:Function}=Plots.plot,
            plotkws::NamedTuple=NamedTuple(), kws...)
        plotkws = (;selectbackground(F)..., plotkws...)
        if fontscale
            kws2 = (;nplots=length(F))
        else
            kws2 = ()
        end
        if as == Dict
            obj = as(key=>show_field(f; kws2..., key=key, kws...) for (key,f) in F)
        elseif as == Plots.plot
            plotkws=(;margin=-2mm, background_color=background,
                     foreground_color=background,
                     background_color_outside=background, plotkws...)
            kws = (;textcolor=textcolor, kws...)
            obj =  as((show_field(f; kws2..., key=key, kws...) for (key,f)
                       in sort(collect(F), by=x->x[1]))...; plotkws...)
        else
            obj =  as(show_field(f; kws2..., key=key, kws...) for (key,f) in F)
        end
        return obj
    end
    function show_field(FF::AbstractVector; func=Plots.bar,
            key::Union{String,NamedTuple,Nothing}=NamedTuple(),
            keyappend::Union{String,NamedTuple,Nothing}=nothing,
            keyprepend::Union{String,NamedTuple,Nothing}=nothing,
            grid=[], textcolor=:black, justification::Symbol=:bottom,
            fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1,
            quant::Vector{Float64}=[0.05,0.99], kws...)
        if all(isnan.(FF)) || isempty(FF)
            return Plots.plot()
        end
        ylims = quantile(utils.skipnan(vec(FF)), quant)
        if grid == []; extrakwargs = (xticks=[], yticks=[])
        else; extrakwargs = ()
        end
        kws = (;extrakwargs..., kws...)
        l = func(grid..., FF; kws...)
        annotate_field(l; key=key, keyappend=keyappend, keyprepend=keyprepend,
                       grid=grid, textcolor=textcolor,
                       justification=justification,
                       location=location,nplots=nplots,
                       nplots_factor=nplots_factor)
    end
    function show_field(FF::AbstractMatrix; 
            key::Union{String,NamedTuple,Nothing}=NamedTuple(),
            keyappend::Union{String,NamedTuple,Nothing}=nothing,
            keyprepend::Union{String,NamedTuple,Nothing}=nothing,
            grid=[], textcolor=:black, justification::Symbol=:bottom,
            fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1,
            quant::Vector{Float64}=[0.05,0.99], kws...)
        if all(isnan.(FF)) || isempty(FF)
            return Plots.plot()
        end
        clims = quantile(utils.skipnan(vec(FF)), quant)
        if grid == []; 
            extrakwargs = (xticks=[], yticks=[])
            grid = (); pos = ();
        else; extrakwargs = (); pos = grid;
        end
        kws = (;extrakwargs..., kws...)
        # ------
        # HEATMAP
        # ------
        hm = heatmap(pos..., FF; clims=Tuple(clims), 
                     colorbar=false, padding=(0,0),
                     c=cgrad(:acton,rev=false), showaxis=:no, kws...)
        #old cms = :linear_kryw_5_100_c67_n256
        #:lajolla
        annotate_field(hm; key=key, keyappend=keyappend, keyprepend=keyprepend,
                       grid=grid, textcolor=textcolor,
                       justification=justification,
                       location=location,nplots=nplots,
                       nplots_factor=nplots_factor)
    end
    function show_field(FF::AbstractArray{T,3}; 
                key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                keyappend::Union{String,NamedTuple,Nothing}=nothing,
                keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                grid=[], textcolor=:black, justification::Symbol=:bottom,
                fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1,
                quant::Vector{Float64}=[0.05,0.99], seriestype=:volume,
                slice_dim=3, kws...) where T <: Real
        if all(isnan.(FF)) || isempty(FF)
            return Plots.plot()
        end
        if seriestype isa Symbol
            @debug "Volume path"
            kws = (; margins=-2mm, padding=(-1,-1), kws...)
                #FF = replace(FF, NaN=>0)
                t=text_annotate(key=key, keyappend=keyappend, 
                                keyprepend=keyprepend,
                               grid=grid, textcolor=textcolor,
                               justification=justification,
                               location=location, nplots=nplots,
                               nplots_factor=nplots_factor)
                kws = (; t..., kws...)
                p = Plots.plot3d(FF; seriestype=seriestype, kws...)
            try
                @debug "volume plotted"
                @debug "annotation plotted"
            catch
                @warn "Failed on key=$key"
            end
        else
            p = mapslices(x->show_field(x;key=key, keyappend=keyappend,
                                    keyprepend=keyprepend, grid=grid,
                                    textcolor=textcolor,
                                    justification=justification,
                                    fontsize=fontsize, location=location,
                                    nplots=nplots,
                                    nplots_factor=nplots_factor,
                                    quant=quant),
                      FF, dims=slice_dim)
        end
        return p
    end
    function annotate_field(hm;
                key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                keyappend::Union{String,NamedTuple,Nothing}=nothing,
                keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                grid=[], textcolor=:black, justification::Symbol=:bottom,
                fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1)
        # ANNOTATION OF KEY
        # -----------------
        fontsize = max(Int(round(fontsize/(nplots*nplots_factor))),3)
        #println("fontsize=$fontsize")
        get_loc(p,l) = p[1] + (p[2]-p[1])*l
        x = get_loc(xlims(), location[1])
        y = get_loc(ylims(), location[2])
        color = textcolor
        font = :bold
        annotation = [transform_key(key)]
        if keyappend != nothing
            push!(annotation, transform_key(keyappend))
        end
        if keyprepend != nothing
            pushfirst!(annotation, transform_key(keyappend))
        end
        annotation = join(annotation, "\n")
        t = text(annotation, justification, color, :right, :bold, pointsize=fontsize)
        annotate!(hm, x, y, t) 
    end
    function text_annotate(;
                key::Union{String,NamedTuple,Nothing}=NamedTuple(),
                keyappend::Union{String,NamedTuple,Nothing}=nothing,
                keyprepend::Union{String,NamedTuple,Nothing}=nothing,
                grid=[], textcolor=:black, justification::Symbol=:bottom,
                fontsize=12, location=(1, 0.03), nplots=1, nplots_factor=1)
        # ANNOTATION OF KEY
        # -----------------
        fontsize = max(Int(round(fontsize/(nplots*nplots_factor))),3)
        #println("fontsize=$fontsize")
        get_loc(p,l) = p[1] + (p[2]-p[1])*l
        x = get_loc(xlims(), location[1])
        y = get_loc(ylims(), location[2])
        color = textcolor
        font = :bold
        annotation = [transform_key(key)]
        if keyappend != nothing
            push!(annotation, transform_key(keyappend))
        end
        if keyprepend != nothing
            pushfirst!(annotation, transform_key(keyappend))
        end
        annotation = join(annotation, "\n")
        return (;title=annotation, titlefontsize=fontsize,
                titlefontcolor=textcolor, titlefontvalign=justification)
    end
    function title_annotate(hm; kws...)
        annotation = text_annotate(;kws...).title
        hm.attr[:plot_title]    = annotation
        hm.attr[:plot_titlefontvaling] = justification
        hm.attr[:plot_titlefontcolor] = color
        hm.attr[:window_title] = annotation
    end
end

export plot
