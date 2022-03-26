
using Blink

function bdisp(item)
    w = Window()
    ui = Interact.@manipulate for i in [1]
        plot(item) 
    end
    body!(w, ui)
end


# ,---.|         |    
# |---'|    ,---.|--- 
# |    |    |   ||    
# `    `---'`---'`---'
# ------------------------------
# General overview: Black Raster
# ------------------------------
raster_overview_black = begin
    set_default_plot_size(20inch, 10inch)
    black_style_vlines = style(point_shapes=[Shape.vline], 
                        default_color=colorant"white",
                        point_size=0.5mm,
                        alphas=[0.5],
                       )
    plot(dspikes, x=:time, y=:unit, black_style_vlines, DefaultGuides...)
end

# ------------------------------
# General overview: simple color
# ------------------------------
colored_by_area = begin 
    set_default_plot_size(20inch, 10inch)
    dspikes.alpha = (dspikes.unit./maximum(dspikes.unit))./2 .+ 0.5;
    plot(dspikes, x=:time, y=:unit, color=:area, alpha=:alpha,
        vertical_lines, DefaultGuides...
       )
end

# ------------------------------
# General overview: rainbow
# ------------------------------

# SET COLORS IN ADVANCE
n(x) = x./maximum(x)
#dspikes[!,"newunit"] = utils.randomize_int(dspikes.unit)
#areacolor = [ColorSchemes.devon.colors, ColorSchemes.lajolla.colors];
areacolor = [ColorSchemes.thermal, ColorSchemes.thermal];

dspikes = groupby(dspikes, "area")
for g = 1:length(dspikes)
    dspikes[g][!,"areaunit2"] = utils.randomize_int(dspikes[g].areaunit)
    if g == 1
        dspikes[g][!, "areaunit_ext"] = dspikes[g].areaunit2
    else
        dspikes[g][!, "areaunit_ext"] = dspikes[g].areaunit2 .+ maximum(dspikes[g-1].areaunit)
    end
    dspikes[g].color_brain_area = get(areacolor[g], n(dspikes[g].areaunit2));
end
dspikes = combine(dspikes, x->x);

# PLOT
set_default_plot_size(20inch, 10inch)
#dspikes.alpha = (dspikes.unit./maximum(dspikes.unit))./2 .+ 0.5;
#dspikes.color_area = ColorSchemes.phase.colors[dspikes.newunit];
#runEpochLayer= layer(epoch, xmin=:start, xmax=:end, color=["gray"], Geom.vband)

getArea(x, area) = x[x.area.==area,:]
ca1 = getArea(dspikes,"CA1")
pfc = getArea(dspikes,"PFC")
#pfcScale = Scale.color_continuous(colormap=x->get(ColorSchemes.diverging_isoluminant_cjm_75_c24_n256, 1-(x/2)))
#ca1Scale = Scale.color_continuous(colormap=x->get(ColorSchemes.diverging_isoluminant_cjm_75_c24_n256, abs.(-x/2)))

for (dark_state, dark_string) in zip((true,false),("dark","light"))

    raster.set_default_theme(dark=dark_state)
    scheme = ColorSchemes.diverging_isoluminant_cjm_75_c24_n256
    scheme = get(scheme, [0, 1])
    scheme = parse.(Colorant, scheme)
    colorManual = Scale.color_discrete_manual(scheme...; levels=["CA1","PFC"])
    neuralLayer = layer(dspikes[begin:1:end,:], vertical_lines, x=:time,
                     color=:area, y=:areaunit_ext, alpha=[0.7])
    coord = Coord.cartesian(xmin=minimum(dspikes.time), xmax=maximum(dspikes.time),
                           ymin=0, ymax=230)
    p = plot(coord, colorManual, neuralLayer, Guide.xlabel("Time"), Guide.ylabel("Neuron"), style(key_position=:none))
    p |> PDF(plotsdir("raster", "raster_ca1pfc_$dark_string.pdf"),dpi=200)
    p |> SVG(plotsdir("raster", "raster_ca1pfc_$dark_string.svg"))
end



# ------------------------------
# General overview: brain area rainbow
# emphasis: instability
# ------------------------------
begin
    dspikes = groupby(dspikes, "area")
    for g = 1:length(dspikes)
        dspikes[g].newareaunit = utils.randomize_int(dspikes[g].areaunit);
        dspikes[g].color_brain_area = areacolor[g][dspikes[g].newareaunit];
    end
    dspikes = combine(dspikes, x->x);
    #dspikes.alpha = (dspikes.unit./maximum(dspikes.unit))./2 .+ 0.5;
    dspikes.newunit = utils.randomize_int(dspikes.unit);
    dspikes.color_area = ColorSchemes.phase.colors[dspikes.newunit];
end

plot(dspikes, x=:time, y=:areaunit, alpha=[0.7],
     ygroup=:area, color=:color_brain_area, vertical_lines,  
     Geom.subplot_grid(Geom.point, free_y_axis=true))

dspikes = begin
    dspikes
end


# -----------------------------
# General overview: brain area rainbow
# emphasis: instability
# ------------------------------
dropmissing!(dspikes)
area_split_instab = plot(dspikes, x=:time, y=:areaunit, alpha=[0.7],
     ygroup=:area, color=:color_brain_area, vertical_lines,  
     Geom.subplot_grid(Geom.point, free_y_axis=true));

# -----------------------------
# General overview: brain area rainbow
# emphasis: stability
# ------------------------------
area_split_stab = plot(dspikes, x=:time, y=:areaunit2, alpha=[0.7],
     ygroup=:area, color=:color_brain_area, vertical_lines,  
     Geom.subplot_grid(Geom.point, free_y_axis=true));


