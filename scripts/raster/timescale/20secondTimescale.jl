
settings["frac_start"] = 0.1
settings["sampling_minutes"] = 0.33
settings["runepoch"] = 3

i = settings["runepoch"];
κ = settings["frac_start"];
epoch = runepoch[i, :];
τ_start = epoch.epoch_start .+ κ*epoch.δ;
τ_end   = τ_start + settings["sampling_minutes"];

timeset = (raster.time .> τ_start) .& (raster.time .< τ_end);
χ = raster[timeset, :];

set_default_plot_size(20inch, 10inch)
χ.newunit = randomize_int(χ.unit)
χ.color_area = ColorSchemes.phase.colors[χ.newunit];
runEpochLayer= layer(runepoch, xmin=:epoch_start, 
                                          xmax=:epoch_end, color=["gray"], Geom.vband)
plot(χ, x=:time, y=:newunit, color=:color_area, alpha=[0.7],
         vertical_lines, DefaultGuides...)

settings["frac_start"] = 0.115
settings["sampling_minutes"] = 0.33
settings["runepoch"] = 3

i = settings["runepoch"];
κ = settings["frac_start"];
epoch = runepoch[i, :];
τ_start = epoch.epoch_start .+ κ*epoch.δ;
τ_end   = τ_start + settings["sampling_minutes"];

timeset = (raster.time .> τ_start) .& (raster.time .< τ_end);
χ = raster[timeset, :];

set_default_plot_size(20inch, 10inch)
χ.newunit = randomize_int(χ.unit)
χ.color_area = ColorSchemes.phase.colors[χ.newunit];
runEpochLayer= layer(runepoch, xmin=:epoch_start, 
                     xmax=:epoch_end, color=["gray"], Geom.vband)
plot(χ, x=:time, y=:newunit, color=:color_area, alpha=[0.7],
    vertical_lines, DefaultGuides...)

