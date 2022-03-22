module goalvector
    using CSV, DataFrames
    using Gadfly
    using Colors, ColorSchemes

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

    """
    TITLE: goalVectorTheme
    Purpose: load gv dataset
    """
    function load()
        X = CSV.read(DrWatson.datadir("exp_pro",
                                      "goal-vector_stops_currentAngle=shuffle.csv"),
                     DataFrame, 
                     strict=false, 
                     missingstring=["NaN", ""]);
        X.stopWell = "Goal " .* string.(X.d_stopWell.-1)
        fr_correction = true # TODO DISSABLE THIS whenever you rerun your tablee extraction
        if fr_correction
            X.occNorm *= 30;
        end
        X = X[filtneuron(X),:]

    end


end
