#+TITLE: Goal-vector Tidy Plots
#+PURPOSE: Takes tidy goal-vector data and visualizes
#+DATE: Jan 24 2022
#
#+-----
#+TODO:
#+-----
#       - Identify list of shuffle significant neurons
#       - Use that list to revisualize the goal-vector tuning (bins × neurons × FR)
#

import DrWatson
DrWatson.quickactivate(expanduser("~/Projects/goal-code/"))

sdps = set_default_plot_size;

# ---------------------
# Angle Related Analysis
# ---------------------

# Main w/o Shuffle
# =================
X = CSV.read(
             DrWatson.datadir("exp_pro",
                              "goal-vector_stops_currentAngle=shuffle.csv"),
             DataFrame, strict=false, missingstring=["NaN", ""])
X.stopWell = "Goal " .* string.(X.d_stopWell.-1)
fr_correction = true # TODO DISSABLE THIS whenever you rerun your tablee extraction
if fr_correction
    X.occNorm *= 30;
end
begin # FR Occ Norm

    # PART 1 : OCCNORM
    # ----------------
    x = X[!, [:d_stopWell, :bins, :d_neuron, :occNorm]]
    uBins = unique(x.bins);
    d = median(diff(uBins)/2);
    eta = 0.1 * d;
    d = d-eta;
    
    # Is there any overall angular distribution?
    # -----
    xr = bootGroups(x, [:bins]);
    goal = plot(xr, xintercept=[0], 
                x=uBins, xmin=uBins.-d, xmax=uBins.+d, y=:mean, 
                ymin = :lower, ymax = :upper,
                Coord.cartesian(xmin=-pi, xmax=pi, ymax=0.4*30),
                Guide.xlabel("Goal Angles"),
                Guide.ylabel("FR Occ Norm"),
                Guide.title("Overall FR Bias"),
                layer(Geom.vline(color="red"), 
                      style(line_style=[:dash], line_width=1mm)),
                layer(Geom.errorbar(), color=[colorant"white"],
                      xmin=:bins, xmax=:bins,
                      style(line_style=[:dash], line_width=1mm)),
                Geom.bar, 
               )

    # Split by goal?
    # -------
    xr = bootGroups(X, [:bins, :d_stopWell, :stopWell])
    layer1 = layer(xmin=xr.bins.-d, xmax=xr.bins.+d, Geom.bar);
    white = repeat([colorant"white"], nrow(xr));
    layer2 = layer(xmin=:bins, xmax=:bins, 
                   Geom.errorbar, 
                   style(line_width=0.5mm, default_color=colorant"white"));
    layer3 = layer(Geom.vline(color="red"), 
                   style(line_width=1mm, line_style=[:dash]));
    subplot = Geom.subplot_grid(layer2, layer3, layer1);
    subplot.coord = Coord.cartesian(xmin=-pi, xmax=pi, ymin=0);
    splitGoal = plot(xr, 
             xgroup=:stopWell, 
             xintercept=[0], x=:bins,
             y=:mean, ymin=:lower, ymax=:upper, 
             Guide.xlabel("Goal Angles"),
             Guide.ylabel("FR Occ Norm"), 
            subplot);


    # Neuron-wise
    # -------
    xr = bootGroups(X, [:bins, :d_neuron, :stopWell])
    groups = groupby(xr, :d_neuron);
    for i in 1:length(groups)
        groups[i].mean = (groups[i].mean .- minimum(groups[i].mean))./(maximum(groups[i].mean)-minimum(groups[i].mean));
    end
    xr = combine(groups, x -> x);
    dropmissing!(xr)
    xr.mean[isnan.(xr.mean)] .= 0
    xr = combine(groupby(xr, :stopWell, sort=true)[2:end], x->x)
    neuronwise = plot(xr, 
             x=:bins, 
             y=:d_neuron, 
             color=:mean,
             xgroup=:stopWell,
             Guide.xlabel("Goal Angles"),
             Guide.ylabel("FR Occ Norm"), 
             Geom.subplot_grid(Geom.rectbin))

    set_default_plot_size(30cm, 20cm)
    angleBiases = vstack(goal, splitGoal)
    set_default_plot_size(25cm, 35cm)
    biasPicture = hstack(angleBiases, neuronwise)

end

begin # FR Occ Norm

    # PART 1 : OCCNORM
    # ----------------
    x = X[filtneuron(X), [:d_stopWell, :bins, :d_neuron, :occNorm]]
    uBins = unique(x.bins);
    d = median(diff(uBins)/2);
    eta = 0.1 * d;
    d = d-eta;
    
    # Is there any overall angular distribution?
    # -----
    xr = bootGroups(x, [:bins]);
    goal = plot(xr, xintercept=[0], 
                x=uBins, xmin=uBins.-d, xmax=uBins.+d, y=:mean, 
                ymin = :lower, ymax = :upper,
                Coord.cartesian(xmin=-pi, xmax=pi, ymax=0.4*30),
                Guide.xlabel("Goal Angles"),
                Guide.ylabel("FR Occ Norm"),
                Guide.title("Overall FR Bias"),
                layer(Geom.vline(color="red"), 
                      style(line_style=[:dash], line_width=1mm)),
                layer(Geom.errorbar(), color=[colorant"white"],
                      xmin=:bins, xmax=:bins,
                      style(line_style=[:dash], line_width=1mm)),
                Geom.bar, 
               );

    # Split by goal?
    # -------
    xr = bootGroups(X[filtneuron(X),:], [:bins, :d_stopWell, :stopWell])
    layer1 = layer(xmin=xr.bins.-d, xmax=xr.bins.+d, Geom.bar);
    white = repeat([colorant"white"], nrow(xr));
    layer2 = layer(xmin=:bins, xmax=:bins, 
                   Geom.errorbar, 
                   style(line_width=0.5mm, default_color=colorant"white"));
    layer3 = layer(Geom.vline(color="red"), 
                   style(line_width=1mm, line_style=[:dash]));
    subplot = Geom.subplot_grid(layer2, layer3, layer1);
    subplot.coord = Coord.cartesian(xmin=-pi, xmax=pi, ymin=0);
    splitGoal = plot(xr, 
             xgroup=:stopWell, 
             xintercept=[0], x=:bins,
             y=:mean, ymin=:lower, ymax=:upper, 
             Guide.xlabel("Goal Angles"),
             Guide.ylabel("FR Occ Norm"), 
            subplot);


    # Neuron-wise
    # -------
    xr = bootGroups(X[filtneuron(X),:], [:bins, :d_neuron, :stopWell])
    groups = groupby(xr, :d_neuron);
    for i in 1:length(groups)
        groups[i].mean = (groups[i].mean .- minimum(groups[i].mean))./(maximum(groups[i].mean)-minimum(groups[i].mean));
    end
    xr = combine(groups, x -> x);
    dropmissing!(xr)
    xr.mean[isnan.(xr.mean)] .= 0
    xr = combine(groupby(xr, :stopWell, sort=true)[2:end], x->x)
    neuronwise = plot(xr, 
             x=:bins, 
             y=:d_neuron, 
             color=:mean,
             xgroup=:stopWell,
             Guide.xlabel("Goal Angles"),
             Guide.ylabel("FR Occ Norm"), 
             Geom.subplot_grid(Geom.rectbin))

    set_default_plot_size(30cm, 20cm)
    angleBiases = vstack(goal, splitGoal)
    set_default_plot_size(25cm, 35cm)
    biasPicture = hstack(angleBiases, neuronwise)

end

# SHUFFLE
# =======
begin
    X = CSV.read(DrWatson.datadir("exp_pro",
                                  "goal-vector_stops_currentAngle=shuffle.csv"),
                 DataFrame, strict=false, missingstring=["NaN", ""])
    X.stopWell = "Goal " .* string.(X.d_stopWell.-1)
    cols = [ "raw_differences", "raw_meanDifferences", "raw_varDifferences", "raw_main", "raw_shuffle", "occNorm_differences", "occNorm_meanDifferences", "occNorm_varDifferences", "occNorm_main", "occNorm_shuffle", "rayleigh_directionality_index_differences", "rayleigh_directionality_index_main", "rayleigh_directionality_index_meanDifferences", "rayleigh_directionality_index_shuffle", "rayleigh_directionality_index_varDifferences", "rayleigh_directionality_index_Z_differences", "rayleigh_directionality_index_Z_main"]
    dropmissing!(X, cols)
    x = X[filtneuron(X), :]
    cols = ["stopWell", "d_stopWell", "bins", "d_neuron", "rayleigh_directionality_index_Z_differences",
        "rayleigh_directionality_index_Z_main", "rayleigh_directionality_index_Z_shuffle"];

    # CURRENT angule distribution
    # --------------------
    begin # distributions
#

        # Density() -Rayleigh Z distros visible overall
        # ---------------------------------
        #+VALUE: HIGH, with neuron filtration
        begin

            set_default_plot_size(30cm,15cm)

            distributions_overall = (
                                     plot(x, alpha=[0.4],
                     x=:rayleigh_directionality_index_Z_main, 
                     Scale.x_log10, Geom.density, 
                        layer(x=:rayleigh_directionality_index_Z_shuffle, 
                              style(default_color=colorant"red"), 
                              Geom.density),
                        Guide.title("All Measurements for each measure M ∈ {Goal × Φ × Cluster}"),
                        Guide.manual_color_key("Measurments",
                                               ["Real", "Shuffle"],
                                               1:2),
                        Guide.xlabel("Rayleigh Ζ"), Guide.ylabel("Density")))

            split_by_goal = (
                plot(x, Scale.x_log10, Scale.y_log10, xgroup=:stopWell,
                        Geom.subplot_grid( layer(x=:rayleigh_directionality_index_Z_main, Geom.density),
                              layer(x=:rayleigh_directionality_index_Z_shuffle, 
                                  style(default_color=colorant"red"), Geom.density)),
                        Guide.title("All Measurements for each measure M ∈ { Φ × Cluster}"),
                        Guide.title("Split by Goal"),
                        Guide.manual_color_key("Measurments",
                                               ["Real", "Shuffle"],
                                               1:2),
                        Guide.xlabel("Rayleigh Ζ"), Guide.ylabel("Density")))

            vstack(distributions_overall, split_by_goal)

        end
        begin # HPC-PFC split

            distributions_overall = (
                                     plot(x, alpha=[0.4],
                     x=:rayleigh_directionality_index_Z_main, 
                     Scale.x_log10, Geom.density,
                        layer(x=:rayleigh_directionality_index_Z_shuffle, 
                              style(default_color=colorant"red"), 
                              Geom.density),
                        Guide.title("All Measurements for each measure M ∈ {Goal × Φ × Cluster}"),
                        Guide.manual_color_key("Measurments",
                                               ["Real", "Shuffle"],
                                               1:2),
                        Guide.xlabel("Rayleigh Ζ"), Guide.ylabel("Density")));

            split_by_goal = (
                plot(x, Scale.x_log10, Scale.y_log10, xgroup=:stopWell,
                        Geom.subplot_grid(
                              layer(x=:rayleigh_directionality_index_Z_main, Geom.density),
                              layer(x=:rayleigh_directionality_index_Z_shuffle, 
                                  style(default_color=colorant"red"), Geom.density)),
                        Guide.title("All Measurements for each measure M ∈ { Φ × Cluster}"),
                        Guide.title("Split by Goal"),
                        Guide.manual_color_key("Measurments",
                                               ["Real", "Shuffle"],
                                               1:2),
                        Guide.xlabel("Rayleigh Ζ"), Guide.ylabel("Density")));

            vstack(distributions_overall, split_by_goal)

        end

        # Histogram() Difference together and split by goal
        # ------------------------------------
        #+VALUE: LOW
        #+COMMENT: These suck. Colored splits better below
        begin
            sdps(20cm, 10cm)
            distributions_overall = (
                     plot(x[filtneuron(x), :], x=:rayleigh_directionality_index_Z_differences, 
                        Geom.histogram(bincount=1000),
                        Coord.cartesian(xmin=-2.5, xmax=3),
                        Guide.title("All Measurements for each measure M ∈ {Goal × Φ × Cluster}"),
                        Guide.xlabel("Rayleigh Ζ"), 
                        Guide.ylabel("Density")
                       )
                   )

            sdps(20cm, 20cm)
            subplot = Geom.subplot_grid(Geom.histogram(bincount=1000, density=true))
            subplot.coord = Coord.cartesian(ymin=0, ymax=0.001)
            split_by_goal = (
                     plot(x[filtneuron(x), :], x=:rayleigh_directionality_index_Z_differences, 
                        xgroup=:stopWell, subplot,
                        Guide.title("All Measurements for each measure M ∈ { Φ × Cluster}"),
                        Guide.title("Split by Goal"),
                        Guide.xlabel("Rayleigh Ζ"), Guide.ylabel("Density")))

        end

        # Visualize percent > shuffle
        # ---------------------------
        #+VALUE: Low
        begin
            # neuron, bin, goal
            x.greater_than = float(x.rayleigh_directionality_index_Z_main .> x.rayleigh_directionality_index_Z_shuffle);
            xr = bootGroups(x, [:bins, :d_neuron]; field=:greater_than)[:, Not(["lower", "upper"])]
            xrn = bootGroups(x, [:d_neuron]; field=:greater_than)[:, Not(["lower", "upper"])]
            xrn.bins .= missing

            xr = unstack(xr,   :bins, :d_neuron, :mean)
            xrn = unstack(xrn, :bins, :d_neuron, :mean)
            xrn = xrn[:, names(xr)]
            sort(names(xr)) == sort(names(xrn))
            (names(xr)) == (names(xrn))
            Matrix((xr .- xrn)[:,2:end])
        end

        # Neurons greater than shuffle
        # and cumulative goal response
        # ---------------------------
        #+VALUE: 7
        #+TODO: [
        #           What causes missing bins ∈ (Goal × Neuron)
        #       ]
        begin

            x = X[:, union(cols[:], ["area"])]
            x.gt_shuffle = x.rayleigh_directionality_index_Z_main .> x.rayleigh_directionality_index_Z_shuffle
            summary = combine(groupby(x, [:stopWell, :bins, :d_neuron]), 
                              :gt_shuffle => (x -> meanboot(x, 500)) => [:mean, :lower, :upper])
            summary.α = zeros(size(summary.lower));
            summary.α .= summary.lower;
            summary.α[summary.α .< 0] .= 0;
            summary.α .= summary.α.^2;
            summary = groupby(summary, :bins)[1]

            # Mean > Shuffle
            # Rectbin
            sNeuron, sGoal = size(summary);
            scale = 2cm
            sdps(0.0025 * sNeuron * scale, 0.5 * sGoal * scale)
            rect_wellSig = plot(summary[filtneuron(summary),:],
                 y=:stopWell, x=:d_neuron, 
                 G[(:x, :stopWell)],
                 G[(:y, :neuron)],
                 Scale.alpha_continuous,
                 alpha=:α,
                 color=:mean, Geom.rectbin)

            # Mean > Shuffle
            # MAKE TWO GOAL-WISE DISTROS, LINEAR AND LOG SCALE
            # Shows concentration to beating the shuffle
            density_goals_color = plot(summary[filtneuron(summary), :],
                    color=:stopWell,
                    G[(:x, :gt_shuffle)],Coord.cartesian(xmin=0, xmax=1),
                    (layer(x=:mean, 
                                            Geom.histogram(bincount=30))));
            density_goals_color_log = plot(summary[filtneuron(summary), :],
                    Scale.y_log2,
                    color=:stopWell,
                    G[(:x, :gt_shuffle)],Coord.cartesian(xmin=0, xmax=1),
                    (layer(x=:mean, 
                                            Geom.histogram(bincount=30))));
            sdps(2*10cm, 18cm)
            Left = hstack([density_goals_color, density_goals_color_log])

            # CUMULATIVE RESPONSES
            unstacked_goal = unstack(summary, 
                                     [:d_neuron, :bins], :stopWell, :mean)
            dropmissing!(unstacked_goal)
            goal_cols = names(unstacked_goal)[end-4:end]
            iter = [0.9, 0.95, 0.99];
            for (i, threshold) in Iterators.enumerate(iter)
                #println(threshold)
                gt = Matrix(unstacked_goal[:, goal_cols] .> threshold)
                goal_portion = convert(Matrix{Int}, gt)
                field_name = string(threshold)
                unstacked_goal[:, field_name] .= sum(goal_portion, dims=2);
            end

            stack_goal = stack(unstacked_goal, Between("0.9","0.99"), 
                               variable_name="Threshold", 
                               value_name="Cumulative_Goals")
            stack_goal[!, "threshold"] = parse.(Float64, stack_goal[!,"Threshold"])
            stack_goal.gt_1 = stack_goal.Cumulative_Goals .>= 1
            stack_goal = groupby(stack_goal, :bins)[1]


            sdps(14cm,20cm)
            coord = Coord.cartesian(ymin=0, ymax=1)
            ticks = Guide.yticks(ticks=LinRange(0, 0.5, 3))
            A = plot(stack_goal, x=:Cumulative_Goals, Scale.x_discrete(levels=0:5),
                 color=:threshold, ygroup=:Threshold,
                 Scale.color_continuous(colormap=Scale.lab_gradient("blue", "red")),
                 Guide.xlabel("Total Goals\nReal>Shuff"),
                 Guide.ylabel("% goal response"),
                 Geom.subplot_grid(coord, ticks, Geom.histogram(density=true)))
            coord = Coord.cartesian(ymin=0, ymax=1)
            ticks = Guide.yticks(label=["0", "0.25", "0.5"],ticks=LinRange(0, 1, 3))
            B = plot(stack_goal, x=:gt_1, Scale.x_discrete(levels=0:1),
                 color=:threshold, ygroup=:Threshold,
                 Scale.color_continuous(colormap=Scale.lab_gradient("blue", "red")),
                 Guide.xlabel("Total Goals\nReal>Shuff"),
                 Guide.ylabel("% goal response"),
                 Geom.subplot_grid(coord, ticks, Geom.histogram(density=true)))
            Right=hstack(A,B)
            
        end
        begin #HPC-PFC

            x = X[:, union(cols[:], ["area"])];
            x.gt_shuffle = x.rayleigh_directionality_index_Z_main .> x.rayleigh_directionality_index_Z_shuffle
            summary = combine(groupby(x, [:stopWell, :bins, :d_neuron, :area]), 
                              :gt_shuffle => (x -> meanboot(x, 500)) => [:mean, :lower, :upper])
            summary.α = zeros(size(summary.lower));
            summary.α .= summary.lower;
            summary.α[summary.α .< 0] .= 0;
            summary.α .= summary.α.^2;

            # Mean > Shuffle
            # Rectbin
            sNeuron, sGoal = size(summary);
            scale = 2cm
            sdps(0.0025 * sNeuron * scale, 0.5 * sGoal * scale)
            rect_wellSig = plot(summary[filtneuron(summary),:],
                 y=:stopWell, x=:d_neuron, 
                 ygroup=:area,
                 G[(:x, :stopWell)],
                 G[(:y, :neuron)],
                 Scale.alpha_continuous,
                 alpha=:α,
                 color=:mean, Geom.subplot_grid(Geom.rectbin))

            # Mean > Shuffle
            # MAKE TWO GOAL-WISE DISTROS, LINEAR AND LOG SCALE
            # Shows concentration to beating the shuffle
            coord = Coord.cartesian(xmin=0, xmax=1);
            density_goals_color = plot(summary[filtneuron(summary), :],
                    color=:stopWell,
                    xgroup=:area,
                    G[(:x, :gt_shuffle)],
                    x=:mean, 
                           Geom.subplot_grid(coord,Geom.histogram(bincount=30)));
            sdps(2*10cm, 18cm)
            density_goals_color

            # CUMULATIVE RESPONSES
            unstacked_goal = unstack(summary, 
                                     [:d_neuron, :bins, :area], :stopWell, :mean)
            dropmissing!(unstacked_goal)
            cols = names(unstacked_goal)[end-4:end]
            iter = [0.9, 0.95, 0.99];
            for (i, threshold) in Iterators.enumerate(iter)
                #println(threshold)
                gt = Matrix(unstacked_goal[:, cols] .> threshold)
                goal_portion = convert(Matrix{Int}, gt)
                field_name = string(threshold)
                unstacked_goal[:, field_name] .= sum(goal_portion, dims=2);
            end

            stack_goal = stack(unstacked_goal, Between("0.9","0.99"), 
                               variable_name="Threshold", 
                               value_name="Cumulative_Goals")
            stack_goal[!, "threshold"] = parse.(Float64, stack_goal[!,"Threshold"])
            stack_goal.gt_1 = stack_goal.Cumulative_Goals .>= 1
            stack_goal = groupby(stack_goal, :bins)[1]


            sdps(14cm,20cm)
            coord = Coord.cartesian(ymin=0, ymax=1)
            ticks = Guide.yticks(ticks=LinRange(0, 0.5, 3))
            A = plot(stack_goal, x=:Cumulative_Goals, Scale.x_discrete(levels=0:5),
                 color=:threshold, ygroup=:Threshold, xgroup=:area,
                 Scale.color_continuous(colormap=Scale.lab_gradient("blue", "red")),
                 Guide.xlabel("Total Goals\nReal>Shuff"),
                 Guide.ylabel("% goal response"),
                 Geom.subplot_grid(coord, ticks, Geom.histogram(density=true)))
            coord = Coord.cartesian(ymin=0, ymax=1)
            ticks = Guide.yticks(label=["0", "0.25", "0.5"],ticks=LinRange(0, 1, 3))
            B = plot(stack_goal, x=:gt_1, Scale.x_discrete(levels=0:1),
                 color=:threshold, ygroup=:Threshold, xgroup=:area,
                 Scale.color_continuous(colormap=Scale.lab_gradient("blue", "red")),
                 Guide.xlabel("Total Goals\nReal>Shuff"),
                 Guide.ylabel("% goal response"),
                 Geom.subplot_grid(coord, ticks, Geom.histogram(density=true)))
            Right=hstack(A,B)
            
        end

        # Visualize well sig
        # ---------------------------
        #+VALUE: 10
        #+COMMENT: Can very clearly see bias above shuffles here and shows
        #          cumulative ceell sig
        begin # ALL CELLS, not split
            not = Not([:rayleigh_directionality_index_Z_main, :rayleigh_directionality_index_Z_shuffle])
            summary = combine(groupby(x[:, not], [:stopWell, :bins, :d_neuron]), 
                    :rayleigh_directionality_index_Z_differences => mean)
            summary = sort(summary, [:stopWell, :d_neuron, :bins])
            summary = groupby(summary, :bins)[1][:, Not(:bins)]
            summary.α = convert(Array{Float32}, (summary.rayleigh_directionality_index_Z_differences_mean .> 1));

            sNeuron, sGoal = size(summary);
            scale = 6cm
            sdps(0.01 * sNeuron * scale, 0.5 * sGoal * scale)
            rect_wellSig = plot(summary,
                 y=:stopWell, x=:d_neuron, 
                 G[(:x, :stopWell)],
                 G[(:y, :neuron)],
                 Scale.alpha_continuous,
                 alpha=:α,
                 color=:rayleigh_directionality_index_Z_differences_mean, Geom.rectbin)


            sdps(8cm, 20cm)
            density_goals_zoom = plot(summary,
                ygroup=:stopWell,
                G[(:x, :rayleighZ_diff)],
                Geom.subplot_grid(layer(x=:rayleigh_directionality_index_Z_differences_mean, 
                                        Geom.histogram(bincount=15))))

            # MAKE TWO GOAL-WISE DISTROS, LINEAR AND LOG SCALE
            # Shows concentration to beating the shuffle
            density_goals_color = plot(summary[filtneuron(summary), :],
                    color=:stopWell,
                    G[(:x, :rayleighZ_diff)],
                    (layer(x=:rayleigh_directionality_index_Z_differences_mean, 
                                            Geom.histogram(bincount=15))));
            density_goals_color_log = plot(summary[filtneuron(summary), :],
                    Scale.y_log2,
                    color=:stopWell,
                    G[(:x, :rayleighZ_diff)],
                    (layer(x=:rayleigh_directionality_index_Z_differences_mean, 
                                            Geom.histogram(bincount=15))));

            sdps(2*8cm, 20cm)
            hstack([density_goals_color, density_goals_color_log])

        end
        begin # Same but HPC-PFC
            x = X[:, union(cols[:], ["area"])]
            not = Not([:rayleigh_directionality_index_Z_main, :rayleigh_directionality_index_Z_shuffle])
            summary = combine(groupby(x[:, not], [:stopWell, :bins, :d_neuron, :area]), 
                    :rayleigh_directionality_index_Z_differences => mean)
            summary = sort(summary, [:stopWell, :d_neuron, :area, :bins])
            summary = groupby(summary, :bins)[1][:, Not(:bins)]
            summary.α = convert(Array{Float32}, (summary.rayleigh_directionality_index_Z_differences_mean .> 1));

            sNeuron, sGoal = size(summary);
            scale = 6cm
            sdps(0.01 * sNeuron * scale, 0.5 * sGoal * scale)
            rect_wellSig = plot(summary,
                 y=:stopWell, x=:d_neuron, 
                 ygroup=:area,
                 G[(:x, :stopWell)],
                 G[(:y, :neuron)],
                 Scale.alpha_continuous,
                 alpha=:α,
                 color=:rayleigh_directionality_index_Z_differences_mean, 
                 Geom.subplot_grid(Geom.rectbin))


            sdps(8cm, 20cm)
            density_goals_zoom = plot(summary,
                ygroup=:stopWell,
                xgroup=:area,
                G[(:x, :rayleighZ_diff)],
                Geom.subplot_grid(layer(x=:rayleigh_directionality_index_Z_differences_mean, 
                                        Geom.histogram(bincount=15))))

            # MAKE TWO GOAL-WISE DISTROS, LINEAR AND LOG SCALE
            # Shows concentration to beating the shuffle
            coord = Coord.cartesian(ymin = 0, ymax=100)
            density_goals_color = plot(summary[filtneuron(summary), :],
                    color=:stopWell,
                    xgroup=:area,
                    x=:rayleigh_directionality_index_Z_differences_mean,
                    G[(:x, :rayleighZ_diff)],
                    (layer(x=:rayleigh_directionality_index_Z_differences_mean, 
                           Geom.subplot_grid(coord, Geom.histogram(bincount=15)))))
            density_goals_color_log = plot(summary[filtneuron(summary), :],
                    Scale.y_log2,
                    color=:stopWell,
                    xgroup=:area,
                    G[(:x, :rayleighZ_diff)],
                    (layer(x=:rayleigh_directionality_index_Z_differences_mean, 
                                            Geom.histogram(bincount=15))));

            sdps(2*8cm, 20cm)
            hstack([density_goals_color, density_goals_color_log])

        end
    end

end


# -------------------------
# Distance related analysis
# -------------------------

# SHUFFLE
# -------
begin
    X = CSV.read(
                 DrWatson.datadir("exp_pro",
                                  "goal-vector_stops_currentDistance=shuffle.csv"),
                 DataFrame, strict=false, missingstring=["NaN", ""])
    X.stopWell = "Goal " .* string.(X.d_stopWell.-1)
    dropmissing!(X)
    cols = ["stopWell", "d_stopWell", "bins", "d_neuron", "rayleigh_directionality_index_Z_differences",
        "rayleigh_directionality_index_Z_main", "rayleigh_directionality_index_Z_shuffle"];
end


# ---------------------
# Saving
# ---------------------
using HDF, JLD
save(plotsdir("goalvector.jld"), "plots", P)
