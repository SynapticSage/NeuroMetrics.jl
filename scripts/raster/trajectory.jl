

settings["frac_start"] = 0.1;
settings["sampling_minutes"] = 2;
settings["runepoch"] = 3;
settings["trajectory_useDownsample"] = false;


# DETERMINE EPOCH PERIOD
i = settings["runepoch"];
κ = settings["frac_start"];
epoch = runepoch[i, :];
τ_start = epoch.start .+ κ*epoch.δ;
τ_end   = τ_start + settings["sampling_minutes"];
# DETERMINE TRAJECTORY TIMES WITHIN
runtraj = get_periods(raster, "traj")
α = findmax(runtraj.end .≥ τ_start)[2];
β = findmin(runtraj.start .≤ τ_end)[2];
trajs = runtraj[ α:β ,:]
trajs = trajs[begin:(end-1),:];
τ_start = minimum(trajs.start);
τ_end = maximum(trajs.end);
# Determine memory periods
memory_periods = _table.binary_on_times(beh, "cuemem", 2,
                                        τ_start=τ_start, 
                                        τ_stop=τ_end)
# Scale pathlength
scale(x, κ) = κ * (x .- minimum(x))./(maximum(x) .- minimum(x));
beh.normCPathLength = scale(beh.currentPathLength, maximum(raster.unit));

# RASTER SNIPPET
if settings["trajectory_useDownsample"] == false
    χ     = _table.constrain_range(raster, "time", τ_start, τ_end)
    χ_beh = _table.constrain_range(beh, "time", τ_start, τ_end)
else
    χ     = _table.constrain_range(downsamp, "time", τ_start, τ_end)
    χ_beh = _table.constrain_range(downsamp_beh, "time", τ_start, τ_end)
end

set_default_plot_size(20inch, 10inch)
χ.newunit = randomize_int(χ.unit)
#χ.color_area = ColorSchemes.phase.colors[χ.newunit];
runTrajLayer = layer(trajs, xmin=:start, 
                     xmax=:end, color=:traj, alpha=[0.35], Geom.vband)
pathLengthLayer = layer(χ_beh, x=:time, y=:x, 
                        color=[colorant"red"], alpha=[0.5],
                        style(line_width=1mm),
                        Geom.line)
vertical_lines = style(point_shapes=[Shape.vline], 
                    line_width=0mm,
                    default_color=colorant"white",
                    point_size=0.5mm)
coordinate = Coord.cartesian(xmin=minimum(χ.time), xmax=τ_end, ymax=maximum(downsamp.unit));
newtrajview = plot(χ, Scale.color_discrete(), runTrajLayer, pathLengthLayer,
     x=:time, y=:newunit, alpha=[0.5], vertical_lines,
     Geom.point, coordinate, DefaultGuides...);
newtrajview
