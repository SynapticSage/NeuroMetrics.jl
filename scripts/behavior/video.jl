using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using GoalFetchAnalysis
import Load, Munge

using Revise, ProgressMeter
using Statistics, NaNStatistics
using StatsPlots: @df
import StatsBase
using VideoIO, GLMakie
using ColorSchemes, Colors
import ColorSchemeTools 
using DataFrames, DataFramesMeta
using Printf, Infiltrator

# --------------------------
# PARAMS
# --------------------------
opt = Dict(
:usevideo                      => false,
:doPrevPast                    => false,
:splitBehVar                   => ["egoVec_1", "egoVec_2", "egoVec_3", "egoVec_4", "egoVec_5"], # Variables to monitor with splits
:vectorToWells                 => true,
:histVectorToWells             => true,
:visualize => :slider,
:nback => 30)

# --------------------------
# LOAD UP THE DATA
# --------------------------
animal, day = "RY16", 36
spikes, beh, tsk  = Load.load(animal, day, data_source=["spikes","behavior","task"])
wells             = Munge.behavior.get_wells_df(animal, day, 2; task=tsk, beh)
boundary          = Munge.task.get_boundary(tsk)
video             = Load.video.load_videocollection(animal, day)
TS = [Load.video.load_videots(animal, day, i) for i in 1:8]
exampframe = video[0.0]
xax, yax = collect.(exampframe.dims)

# --------------------------
# FILTER BY SOME CONDITIONS?
#[ --------------------------
beh = @subset(beh, :epoch .== 2)
T = size(beh, 1)
using CategoricalArrays

function barvar(Fig, b, sym; i = size(Fig.layout)[2])
    bsym = @lift ($b)[sym]
    stringtype = String <: eltype(typeof(beh[:, sym]))
x,y = nothing, nothing
    if !stringtype
        x =  Float64.(unique(skipmissing(beh[:,sym])))
    else
        x = sort(CategoricalArray(unique(skipmissing(beh[:,sym]))))
    end
    ax = stringtype ? 
    Axis(Fig[i,3], xticks=levelcode.(x), yticks=[0,1], yticklabels=("off","on"), title=string(sym), 
         xticklabels=levels(x)) :
    Axis(Fig[i,3], xticks=x, yticks=[0,1], yticklabels=("off","on"), title=string(sym), 
         )
    @infiltrate stringtype
    y = @lift ismissing($bsym) ? zeros(size(x)) : Float64.($bsym .== x)
    GLMakie.current_axis!(ax)
    GLMakie.barplot!(stringtype ? levelcode.(x) : x, y)
end

begin
    N = 5
    Fig = Figure(resolution=(800,800))
    if opt[:visualize] == :video
        t = Observable(Int(1000))
    elseif opt[:visualize] == :slider
        sl_t = Slider(Fig[1:N, 4], range = 1:size(beh,1), horizontal = false, 
                      startvalue = 1000)
        t = lift(sl_t.value) do tt
            tt
        end
    end
    vidax = Axis(Fig[1:N, 1:2])
    B = @lift beh[max(1,$t-opt[:nback]):$t,:]
    b = @lift beh[$t, :]
    # Track and Video
    begin
        frame = @lift Matrix(video[$b.time])'
        realvideo = image!(vidax, yax,xax, frame, zindex=-10, alpha=0.5)
        GLMakie.lines!(vidax, eachcol(boundary)...; color=:red)
    end
    #Behavior
    begin
        Bx, By = @lift($B.x),@lift($B.y)
        bx, by = @lift([$b.x]),@lift([$b.y])
        GLMakie.scatter!(vidax, Bx, By)
    end
    # Reporters
    begin
        barvar(Fig, b, :startWell; i=1)
        barvar(Fig, b, :stopWell;  i=2)
        barvar(Fig, b, :cuemem;    i=3)
        barvar(Fig, b, :correct;   i=4)
        barvar(Fig, b, :hatraj;    i=5)
    end
    #GLMakie.scatter!(vidax, bx, by, c=:orange)
    #cuemem, corr, hatraj, hatrajnum
    #well,   poke
    Fig
end

for o in t[]:(T-100)
    t[] = Int(o)
    sleep(1/90)
end
