using Statistics
using Plots
using Measures
using Printf
using ProgressMeter
quickactivate(expanduser("~/Projects/goal-code"))
includet(srcdir("raw.jl"))
includet(srcdir("behavior.jl"))

day, epoch = 9, 2
paths = raw.load_pathtable("RY22",day)
position = raw.dlc.load("RY22", day, epoch)
crop_x, crop_y = extrema(sub.X), extrema(sub.Y)
usevideo = true
if usevideo
    video = raw.video.load("RY22", day, epoch)
end
path = behavior.path

function name(type, pathsubset)
    name = "day=$(@sprintf("%02d",day)),epoch=$(@sprintf("%02d",epoch)),block=$(pathsubset.block[1])"
    folder = plotsdir("behavior", "dynamic", type)
    if !(isdir(folder))
        mkdir(folder)
    end
    name = plotsdir(folder, name)
end


uBlocks = unique(paths.block)
@showprogress for block in uBlocks
    pathsubset = paths[paths.block.==block, :]
    println(block, " ", size(pathsubset))
    overall, indiv = path.subplot_block(position, pathsubset, video=video); p[1]
    n = name("(B) subplot_per_block_of_traj", pathsubset)
    savefig(overall, n*".png")
    savefig(overall, n*".pdf")
    savefig(overall, n*".svg")
    for i in 1:length(indiv)
        n = name("individual_trajectories", pathsubset)
        savefig(indiv[i], n*"_traj=$(@sprintf("%02d",i)).png")
        savefig(indiv[i], n*"_traj=$(@sprintf("%02d",i)).pdf")
        savefig(indiv[i], n*"_traj=$(@sprintf("%02d",i)).svg")
    end
    p = path.plot_block(position,    pathsubset, video=video); p

    n = name("(A) overlaid_trajectory_blocks", pathsubset)
    savefig(p, n*".png")
    savefig(p, n*".pdf")
    savefig(p, n*".svg")
end
