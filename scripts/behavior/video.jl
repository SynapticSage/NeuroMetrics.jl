using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using GoalFetchAnalysis
import Load, Munge

using Revise, ProgressMeter, CategoricalArrays, Statistics, NaNStatistics, VideoIO, GLMakie, ColorSchemes, Colors, DataFrames, DataFramesMeta, Printf, Infiltrator, ImageFiltering
using StatsPlots: @df
import StatsBase, ColorSchemeTools 

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
:nback => 10)

# --------------------------
# LOAD UP THE DATA
# --------------------------
animal, day = "RY16", 36
#animal, day = "RY22", 21
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


# FUnctions
begin
    function mymktemp()
        path = mktemp()[1] * ".mp4"
        expanduser(replace(path, "/tmp/" => "~/tmp/"))
    end
    function barvar(Fig, b, sym; i = size(Fig.layout)[2])
        bsym = @lift ($b)[sym]
        stringtype = String <: eltype(typeof(beh[:, sym]))
    x,y = nothing, nothing
        if !stringtype
            x =  Float64.(unique(skipmissing(beh[:,sym])))
            y = @lift ismissing($bsym) ? zeros(size(x)) : Float64.($bsym .== x)
        else
            x = sort(CategoricalArray(unique(skipmissing(beh[:,sym]))))
            y = @lift ismissing($bsym) ? zeros(size(x)) : Float64.($bsym .== levels(x))
        end
        ax = Axis(Fig[i,3], title=string(sym))
        if stringtype
            ax.xticks=(levelcode.(x), levels(x))
            ax.yticks=(0:1, ["off","on"])
        else
            ax.xticks=x
            ax.yticks=(0:1, ["off","on"])
        end
        @infiltrate stringtype
        GLMakie.current_axis!(ax)
        GLMakie.barplot!(stringtype ? levelcode.(x) : x, y)
    end
end

# Setup observables
begin

    N = 4
    Fig = Figure(resolution=(1400,700))
    if opt[:visualize] == :video
    t = Observable(Int(1000))
    elseif opt[:visualize] == :slider
    sl_t = Slider(Fig[1:N, 4], 
                  range = 1:size(beh,1), 
                  horizontal=false, 
                  startvalue=1000)
    t = lift(sl_t.value) do tt;
        tt; end
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
        movement = @lift $b.moving
        Bx, By, L = @lift($B.x),@lift($B.y), 
                        @lift collect(LinRange(0,1,length($B.y)))
        bx, by = @lift([$movement ? $b.x : NaN]),@lift([$movement ? $b.y : NaN])
        GLMakie.scatter!(vidax, Bx, By)
        GLMakie.scatter!(vidax, bx, by, color=RGBA(1,0.5,1,0.15), markersize=80)
        hatraj = @lift $b.hatraj === missing ? [""] : [uppercase($b.hatraj)]
        ymid = mean(vidax.yaxis.attributes[:limits][])
        xmid = mean(vidax.xaxis.attributes[:limits][])
        #a=GLMakie.annotations!(hatraj, [Point2f((135.0, 50.0))], color=:orange, fontsize=40)
        cuemem_color = @lift if $b.cuemem == 1
            :orange
        elseif $b.cuemem == 0
            :blue
        else
            :white
        end
        a=GLMakie.annotations!(hatraj, [Point2f((xmid, ymid))], color=cuemem_color, fontsize=40)
        animal_location = @lift [Point2f($b.x,$b.y)]
        speed = @lift round(abs($b.velVec),sigdigits=2)
        speedtext = @lift ["| v⃗ | $($speed)"]
        vtext=GLMakie.annotations!(speedtext, animal_location, color=:hotpink, fontsize=20, align=(:left,:top))
        smspeed = @lift round(abs($b.smoothvel),sigdigits=2)
        smspeedtext = @lift ["sm( v⃗ ) $($smspeed)"]
        smvtext=GLMakie.annotations!(smspeedtext, animal_location, color=:pink, fontsize=20, align=(:left,:top), offset=(0,30))
    end

    begin
        poke = @lift $b.poke
        wellpoke = @lift $poke == 0 ? (NaN,NaN) : Tuple(wells[$poke, [:x,:y]])
        GLMakie.scatter!(vidax, wellpoke, color=:orange, glowwidth=5, glowcolor=:orange, transparency=true)
    end

    # Reporters
    begin
        barvar(Fig, b, :startWell; i=1)
        barvar(Fig, b, :stopWell;  i=2)
        barvar(Fig, b, :cuemem;    i=3)
        barvar(Fig, b, :correct;   i=4)
        #barvar(Fig, b, :hatrajnum; i=5)
    end
    # Well scatter

end
display(Fig)

#@showprogress for stamp in t[]:(T-100)
#    try
#    	t[] = stamp
#    catch
#    	@warn "timestamp failed with t=$(t[])"
#    end
#    sleep(1/90)
#end

prog = Progress(length(t[]:(T-100)), desc="Recording video")
iter = ((next!(prog); t) for t in t[]:(T-100))
GLMakie.record(Fig, mymktemp(), framerate=60, iter, compression=12) do stamp
    try
            t[] = stamp
    catch
            @warn "timestamp failed with t=$(t[])"
    end
end
