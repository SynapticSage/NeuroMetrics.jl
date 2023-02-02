using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using GoalFetchAnalysis, Munge.causal, Munge.manifold
import Load, Munge

using Revise, ProgressMeter, CategoricalArrays, Statistics, NaNStatistics, VideoIO, GLMakie, ColorSchemes, Colors, DataFrames, DataFramesMeta, Printf, Infiltrator, ImageFiltering, TextWrap
using StatsPlots: @df
import StatsBase, ColorSchemeTools 


## ----------
## PARAMETERS
## ----------
animal, day, filt, N = "RY16", 36, :all, 5
areas = (:ca1,:pfc)
distance = :many
feature_engineer = :many # many | nothing

opt = Dict(
:usevideo                      => false,
:doPrevPast                    => false,
:visualize => :slider,
:nback => 10)


# Load all requisite vars
manifold.load_manis_workspace(Main, animal, day; filt, 
      areas, distance, feature_engineer, 
      N)

spikes, beh, ripples, cells  = Load.load(animal, day)
storage = load_alltimes_savefile(animal, day, N; params)
tsk = Load.load_task(animal, day)
wells             = Munge.behavior.get_wells_df(animal, day, 2; task=tsk, beh)
boundary          = Munge.task.get_boundary(tsk)
video             = Load.video.load_videocollection(animal, day)
TS = [Load.video.load_videots(animal, day, i) for i in 1:8]
exampframe = video[0.0]
xax, yax = collect.(exampframe.dims)
if (:index ∉ propertynames(beh))
    beh[!,:index] = 1:size(beh,1)
end

# -------------------------
# Filter by some conditions
# -------------------------
# Behavior
beh = @subset(beh, :epoch .== 2)
T = size(beh, 1)

# Filter manifold
emdfs = @subset(emdf, 
              :metric  .== Symbol("Euclidean"),
              :feature .== Symbol("zscore"),
              :dim .== 3)
isempty(emdfs ) ? @error("emdfs is empty!") : nothing

# Throw away unused partitions
begin
end
# Create a sorted set of partitions
begin
	keynames = [:min_dist, :n_neighbors, :metric, :area, :dim, :feature]
        emdfs = groupby(emdfs, keynames)
end
isempty(emdfs ) ? @error("emdfs is empty!") : nothing
beh.T_partition = Utils.searchsortedprevious.([sort(unique(emdfs.T_start))], beh.index)

println("Number of embedding keys ", size(emdfs,1))
#println("Estimated number of manifolds to plot ", size(emdfs,1)/maximum(emdfs[!,:T_partition]))
println("Estimated number of manifolds to plot ", length(emdfs))

# Align each manifold?


# ----------------
# Helper Functions
# ----------------
begin
    function mymktemp()
        path = mktemp()[1] * ".mp4"
        expanduser(replace(path, "/tmp/" => "~/tmp/"))
    end
end

# -----------------
# Setup observables
# -----------------
begin

    N = length(emdfs)
    Fig = Figure(resolution=(1400,700))
    if opt[:visualize] == :video
    t = Observable(Int(1000))
    elseif opt[:visualize] == :slider
    sl_t = Slider(Fig[1:N, 4], 
                  range = 1:size(beh,1), 
                  horizontal=false, 
                  startvalue=1000)
    sl_phi = Slider(Fig[1:N, 4], 
                  range = 1:size(beh,1), 
                  horizontal=false, 
                  startvalue=1000)
    t = lift(sl_t.value) do tt;
        tt; end
    end
    # can I map this back

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
        poke = @lift $b.poke
        wellpoke = @lift $poke == 0 ? (NaN,NaN) : Tuple(wells[$poke, [:x,:y]])
        GLMakie.scatter!(vidax, wellpoke, color=:orange, glowwidth=5, glowcolor=:orange, transparency=true)
    end

    # Instatiate current overall manifold
    m,M = Dict(),Dict()
    partition = Dict()
    maxes = Dict()

    begin

        # Plot a light grey of the total manis
        for (i, manis) in enumerate(emdfs)
                title = NamedTuple(col=>value for (col,value) 
                                                   in zip(emdfs.cols, inv(emdfs.keymap)[i]))
                title = TextWrap.wrap(replace(Utils.namedtup.tostring(title),","=>" "),
                                                          width=40)
                maxes[i] = Axis3(Fig[i,3], title=title, titlesize=8)
                part = @lift findfirst($t .>= manis.T_start .&& $t .< manis.T_end)		
                partition[i] = @lift $part === nothing ? @error("partition missing") : $part
                e = []
                for dim in range(1, manis.dim[1])
                        push!(e, @lift manis.value[$(partition[i])][:,dim])
                end
                black=RGBA(colorant"black", 0.01)
                GLMakie.scatter!(e..., color=black, transparency=true, markersize=4)
        end

    end

    # Record current manifold point
    begin

        # Plot current location of manis
        for (i, manis) in enumerate(emdfs)
                mani     =  @lift manis[$(partition[i]), :]
                mani_ind = @lift findfirst($mani.inds_of_t .== $b.index)
                mani_start = @lift max($mani_ind-opt[:nback], $mani.T_end)
                point = @lift Point3f($mani.value[$mani_ind,:])
                points = @lift Point3f.(eachrow($mani.value[$mani_start:$mani_ind,:]))
                s=GLMakie.scatter!(maxes[i], point, c=colorant"red",
                                   markersize=6, glowwidth=5,
                                   glowcolor=colorant"red", )
                S=GLMakie.scatter!(maxes[i], point, c=colorant"red",
                                   markersize=6, glowwidth=5,
                                   glowcolor=colorant"red", )
        end
    end

end
display(Fig)

@showprogress for stamp in t[]:(T-100)
    try
    	t[] = stamp
    catch
    	@warn "timestamp failed with t=$(t[])"
		sleep(0.25)
    end
    sleep(1/90)
end

prog = Progress(length(t[]:(T-100)), desc="Recording video")
iter = ((next!(prog); t) for t in t[]:(T-100))
GLMakie.record(Fig, mymktemp(), framerate=60, iter, compression=12) do stamp
    try
            t[] = stamp
    catch
            @warn "timestamp failed with t=$(t[])"
    end
end
