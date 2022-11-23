"""
`rf`
"""
module receptivefield

    using Plots, LaTeXStrings, Measures
    using DataFramesMeta
    using Statistics
    using ProgressMeter
    using NaNStatistics
    import TextWrap
    using Infiltrator
    using Interpolations
    using ImageFiltering
    using Measures
    using Colors

    import ..Plot
    using Field
    using Field: ReceptiveField, ReceptiveFields
    using Utils.binning
    import Utils

    import Random
    using Memoization
    function clear_memoize_receptivefield!() 
        Memoization.empty_cache!(upsample)
        Memoization.empty_cache!(gaussianfilt)
    end
    export clear_memoize_receptivefield!

    struct BarDisplay
        x::Float64
        y::Float64
        xw::Float64
        yw::Float64
        fraction::Float64
        text::String
        foreground::Colorant
        background::Colorant
    end

    """
    Statistic filled bar that appears inset within a plot
    """
    @recipe function plot_bardisplay(bd::BarDisplay)
    end

    @userplot StatBarAx
    function plot(plt::StatBarAx)
        Y = ylims()
        X = xlims()
    end

    @recipe function plot_adaptiverf(field::ReceptiveField, val::Symbol=:rate,
            transpose::Bool=true;
            ztransform::Bool=false, mfunc::Function=nanmean,
            sfunc::Function=nanstd, title_width=40, 
            upsamp::Int=2, gauss::Real=nothing, 
            bardisp=nothing,
            interpmode=(BSpline(Constant()), BSpline(Constant())))
        zz = copy(getproperty(field, val))
        if ndims(zz) == 2
            seriestype --> :heatmap 
        elseif ndims(zz) == 1
            seriestype --> :line
        elseif ndims(zz) == 3
            seriestype --> :volume
        end

        title --> string(field.metrics; width=title_width)
        xx = [field.grid.centers[1]...]
        x --> xx
        if ztransform
            zz = (zz .- mfunc(zz))./sfunc(zz)
            colorbar_title --> String(val) * " Z-transform"
        else
            colorbar_title --> String(val)
        end

        if gauss != 0 && gauss != (0,0) ||
            (gauss === nothing && upsamp != 1)
            zz = gaussianfilt(zz; gauss, upsamp)
        end

        if length(field.grid.centers) > 1
            aspect_ratio --> 1
            yy = [field.grid.centers[2]...]
            y --> yy
            if upsamp != 1
                xx,yy,zz = upsample(xx,yy,zz; interpmode, upsamp)
            end
            transpose ? (xx, yy, zz') : _transpose(xx, yy, zz')
        elseif length(field.grid.centers) == 1
            @series begin
                seriestype --> :line
                label --> "firing " *string(val)
                #fillrange --> (0,zz)
                xlim --> (minimum(xx), maximum(xx))
                xx, zz
            end
            @series begin
                seriestype := :vline
                c := :black
                linestyle := :dash
                label := ""
                ([0])
            end
            @series begin
                seriestype := :scatter
                markersize := 4
                label := ""
                c := :black
                #fillrange := nothing
                xx, zz
            end
            xlim := (nanminimum(xx), nanmaximum(xx))
            ylim := (nanminimum(zz)*0.95, nanmaximum(zz)*1.1)
            label := ""
            ([NaN],[NaN])
        end
    end

    @memoize function gaussianfilt(zz; gauss, upsamp)
        gauss = gauss === nothing ? upsamp * 0.33 : gauss
        gauss = gauss isa Tuple ? gauss : Tuple(gauss for _ in 1:ndims(zz))
        nanzz = isnan.(zz)
        zz = replace(zz, NaN=>0)
        zz = imfilter(zz, Kernel.gaussian(gauss))
        zz[nanzz] .= NaN
        zz = reshape(zz, size(nanzz))
    end

    @memoize function upsample(xx, yy, zz; interpmode, upsamp)
        if interpmode isa Function || interpmode isa UnionAll ||
           interpmode isa Type
            interpmode = (interpmode(Linear()), interpmode(Linear()))
        end
        int = interpolate(zz, interpmode)
        axs = (1:(1/upsamp):size(zz,1), 1:(1/upsamp):size(zz,2))
        zz = int[axs...]
        int = interpolate(xx, BSpline(Linear()))
        axs = 1:(1/upsamp):length(xx)
        xx = int[axs]
        int = interpolate(yy, BSpline(Linear()))
        axs = 1:(1/upsamp):length(yy)
        yy = int[axs]
        xx, yy, zz
    end

    @recipe function plot_adaptiveocc(field::T where T<:Occupancy,
            val::Symbol=:prob, transpose::Bool=true)
        @info "tet"
        a=1
        Z = if val==:prob
            reshape(getproperty(field, val), size(field.grid))
        else
            getproperty(grid, val)
        end
        z := Z
        if ndims(Z) == 2
            seriestype --> :heatmap 
        elseif ndims(Z) == 1
            seriestype --> :line
        elseif ndims(Z) == 3
            seriestype --> :volume
        end
        colorbar_title --> String(val)
        X = [field.grid.centers[1]...]
        x --> X
        if length(field.grid.centers) > 1
            Y = [field.grid.centers[2]...]
            y --> Y
            transpose ? (X,Y,Z') : _transpose(X,Y,Z')
        else
            label --> "occupancy " * string(val)
            fillrange --> (0,Z)
            ylim --> (0, maximum(Z)*1.3)
            (X,Z)
        end
    end

    @recipe function plot_adaptivegrid(grid::GridAdaptive, val::Symbol=:radii,
        transpose::Bool=true)

        Z = getproperty(grid, val)
        Z = Z[1] isa AbstractVector ? sqrt.(prod.(Z)) : Z

        colorbar_title --> String(val)
        if ndims(Z) == 2
            seriestype --> :heatmap 
        elseif ndims(Z) == 1
            seriestype --> :line
        elseif ndims(Z) == 3
            seriestype --> :volume
        end
        c --> :thermal
        X = [grid.centers[1]...]
        x --> X
        if length(grid.centers) > 1
            Y = [grid.centers[2]...]
            y --> Y
            transpose ? (X,Y,Z') : _transpose(X,Y,Z')
        else
            label --> string(val)
            fillrange --> (0,Z)
            ylim --> (0, maximum(Z)*1.3)
            (X,Z)
        end
    end

    function _transpose(X, Y, Z)
        (Y,X,Z')
    end

    
    @userplot FieldSplay
    @recipe function fieldsplay(plt::FieldSplay, val::Symbol=:rate; 
            srt=nothing, vectorize=true, veccolwidth=nothing,
            area_conditional_c::Union{Nothing, AbstractDict}=nothing,
            scale_spg=1, randsamp=nothing, totalsamp=nothing,
            change_on_metric=nothing,
            transpose=true,
            sizepergraph=(150,150), totalarea_aspect=1, metricfilter=nothing,
            metricdisp=nothing
        )
        fields = plt.args[1]
        @assert fields isa ReceptiveFields
        fields = deepcopy(fields)
        if metricfilter isa Function
            fields = filter(x->metricfilter(x.metrics), fields)
        end
        
        if srt !== nothing
            srt = !(srt isa Vector) ? [srt] : srt
            metvals = [field.metrics[met] for field in fields,
                       met in srt]
            metvals = [metvals 1:size(metvals,1)]
            metvals = sortslices(metvals, dims=1,
                       by=x->Tuple(x[i] for i in 1:size(metvals,2)), rev=true)
            fields = fields[Int.(metvals[:,end])]
        end

        fields = totalsamp === nothing ? fields : fields[1:totalsamp]
        fields = randsamp === nothing ? fields : Random.shuffle(fields)[1:randsamp]


        blanks, usefieldsize, undefined = [], true, 0
        veccolwidth = veccolwidth === nothing ?
        round(ceil(sqrt(length(fields))) * totalarea_aspect) :
            veccolwidth
        @assert veccolwidth !== nothing
        if vectorize
            fields = vec(fields)
            undefined = veccolwidth-mod(length(fields), veccolwidth)
            blanks = [nothing for _ in 1:undefined]
            usefieldsize=false
        end
        if metricdisp !== nothing
            for (f,field) in enumerate(fields)
                [pop!(field.metrics, key) for key in keys(field.metrics)
                    if key ∉ metricdisp]
            end
        end
        #fields = [plot(field, val) for field in fields]
        fields = isempty(blanks) ? fields : [fields; blanks]
        gridsize = usefieldsize ? (size(fields)...,) :
            (Int(round(ceil(length(fields))/veccolwidth)),Int(veccolwidth),)
        layout := Plots.grid(gridsize...)

        finalsize = scale_spg .* gridsize .* sizepergraph
        @info "sizes" gridsize size(fields) finalsize
        size --> finalsize
        title_font_pointsize --> 7/finalsize[1]
        #background_color_outside := ;match
        aspect_ratio --> 1


        for f in 1:length(fields)
            @series begin
                margin --> 0.2mm
                subplot --> f
                legend --> false
                grid --> false
                framestyle --> :none
                if length(fields)-f < undefined
                    background_color_inside := :match
                    legend := false
                    grid := false
                    framestyle := :none
                    ()
                else
                    if change_on_metric !== nothing
                        for (change_detect,change_dict) in change_on_metric
                            if change_detect(fields[f])
                                change_dict(d)
                            end
                        end
                    end
                    if area_conditional_c !== nothing
                        c --> area_conditional_c[fields[f].metrics[:area]]
                    end
                    (fields[f], val, transpose)
                end
            end
        end
        ()
    end

end
