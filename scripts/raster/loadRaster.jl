
    # And remove the minimum time, to count from 0
    # Divide, converting to minutes
#+PURPOSE: Raster plot of neural data for grant
#+DATE: 8-2-22
#
# ----------------------------------------------
# TODO
# ----------------------------------------------
# > By area unit labels
# > Randomize unit labels
# > EPOCH hspan geometry

using DrWatson
quickactivate(expanduser("~/Projects/goal-code/"))
cd(DrWatson.projectdir())
using Revise
includet(srcdir("table.jl"))
includet(srcdir("raw.jl"))
includet(srcdir("raster.jl"))
includet(srcdir("utils.jl"))
using Gadfly, Colors, ColorSchemes
using CSV, DataFrames
using Shuffle

#                                                      
#,---.               o                  |     |         
#|---|,---.,---..   ..,---.,---.    ,---|,---.|--- ,---.
#|   ||    |   ||   |||    |---'    |   |,---||    ,---|
#`   '`---'`---|`---'``    `---'    `---'`---^`---'`---^
#              |                                        
# SETTINGS
settings = Dict(
            "downsample"=>true,
            "dfactor"  => 10,
           )
spikes, behavior  = raw.load("RY16", 36);
if settings["downsample"]
    dspikes      = raw.downsample(spikes, dfactor=settings["dfactor"]);
    dbeh = raw.downsample(behavior, dfactor=settings["dfactor"]);
end

#                                                                       
# ,--.      ,---.          |    |        o                              
# |   |,---.|__. ,---..   .|    |---     .,-.-.,---.,---.,---.,---.,   .
# |   ||---'|    ,---||   ||    |        || | |,---||   ||---'|    |   |
# `--' `---'`    `---^`---'`---'`---'    `` ' '`---^`---|`---'`    `---|
#                                                   `---'          `---'
# THEME
raster.set_default_theme();

# ELEMENTS
DefaultGuides = [Guide.xlabel("Time\n(Minutes)"), 
                 Guide.ylabel("Unit"),
                 Guide.title("Recording Day 36, RY16"),
                 Coord.cartesian(ymax=maximum(dspikes.unit))
                ]
DefaultGuidesSubplot = [Guide.xlabel("Time\n(Minutes)"), 
                 Guide.ylabel("Unit"),
                 Guide.title("Recording Day 36, RY16"),
                ]
vscale = 0.1mm
vertical_lines = style(point_shapes=[Shape.vline], 
                    highlight_width=0.1 * vscale,
                    line_width=0mm,
                    point_size=vscale)
areacolor = [ColorSchemes.devon.colors, ColorSchemes.lajolla.colors];

epoch = table.get_periods(spikes, "epoch")
table.remove_interperiod_time!(dspikes, epoch)


