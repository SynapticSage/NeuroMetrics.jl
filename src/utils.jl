module utils

import Random
using CSV, DataFrames
using Gadfly
using Colors, ColorSchemes

export skipnan
export itsizeof, piso

skipnan(x) = Iterators.filter(!isnan, x)

function itsizeof(X)
    [size(x) for x in X]
end
function piso(X)
    println([size(x) for x in X])
end

function randomize_int(X)
    Xmin = minimum(X);
    Xmax = maximum(X);
    initial = collect(Xmin:Xmax);
    final   = Random.shuffle(initial);
    mapping(x) = Dict(initial .=> final)[x]
    map(mapping, X)
end

guides = Dict(); # shortcut for guides
guides[(:x, :stopWell)] = Guide.xlabel("Goal")
guides[(:y, :stopWell)] = Guide.ylabel("Goal")
guides[(:x, :neuron)] = Guide.xlabel("Neuron")
guides[(:y, :neuron)] = Guide.ylabel("Neuron")
guides[(:x, :rayleighZ)] = Guide.xlabel("Rayleigh Ζ")
guides[(:x, :rayleighZ_diff)] = Guide.xlabel("Rayleigh Z\nDifferences")
guides[(:x, :gt_shuffle)] = Guide.xlabel("Percent\nreal > shuff") 

"""
TITLE: goalVectorTheme
Purpose: theme for goal vector shit
"""
function goalVectorTheme()
    theme = Theme(major_label_color=colorant"white", major_label_font_size=14pt,
                  minor_label_color=colorant"white",
                  key_label_color=colorant"white",
                  key_title_color=colorant"white",
                  panel_fill=colorant"black",
                  background_color=colorant"black")
    Gadfly.push_theme(theme)
end

end
