using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
include(scriptsdir("decode","Initialize.jl"))
load_from_checkpoint = true
dothresh = false
dosweep = false
doPrevPast = false
doThetaPhase, doRipplePhase = false, false
animal,day,epoch="RY16",36 ,2

if !(load_from_checkpoint)

    include(scriptsdir("decode", "LoadData.jl"))

    # --------------------
    # Preprocess BEHAHVIOR
    # TODO Turn this into a function
    # --------------------
    beh = Munge.behavior.annotate_pastFutureGoals(beh; doPrevPast)

    # --------------
    # Preprocess LFP
    # --------------
    include(scriptsdir("decode", "PreprocessLFP.jl"))

    Munge.spiking.isolated(spikes, lfp)

    # Checkpoint pre-video data
    Decode.save_checkpoint(Main, decode_file; split=0)
else

    D = Decode.load_checkpoint("/Volumes/FastData/decode/sortedspike.empirical.notbinned.n_split=4.downsamp=1.speedup=20.0/split=0_decode.h5")
    for (key,value) in D
        eval(Meta.parse("$key = D[:$key]"))
    end
    for df in [beh, lfp, cycles, ripples]
        Load.fix_complex(df)
        Load.fix_rgba(df)
    end

end


# -------- END OF CHECKPOINTED DATA PREPROC ----------------------------

# GET A DECODE OBJECT PER WAVE BAND OF INTEREST
lfp, theta, ripples, non = begin
    # ------------------------------
    # Separate DECODES by EVENTS
    # ------------------------------
    using Munge.lfp_decode
    @time theta, ripple, non = separate_theta_ripple_and_non_decodes(T, lfp, dat;
                                                               quantile_range=(0.9,1),
                                                               doThetaPhase,
                                                               doRipplePhase);
    if dosweep
        theta, ripples = convert_to_sweeps(lfp, theta, ripples; doRipplePhase)
    end
    no_cycle_happening = Utils.squeeze(all(isnan.(theta), dims=(1,2)))
    lfp[!, :color_phase] = RGBA.(get(ColorSchemes.romaO, lfp[!, :phase]))
    lfp[Utils.searchsortednearest.([lfp.time], T[no_cycle_happening]), :color_phase] .= RGBA{Float64}(0, 0, 0, 0)
    lfp, theta, ripples, non
end

# RELOAD WELL INFO
@time task   = Load.load_task(animal,     day)
wells = task[(task.name.=="welllocs") .& (task.epoch .== epoch), :]
homeWell, arenaWells = begin
    hw = argmax(StatsBase.fit(StatsBase.Histogram, filter(b->b!=-1,beh.stopWell), 1:6).weights)
    @info "Homewell = $hw"
    wells[hw,:], wells[setdiff(1:5,hw),:]
end
boundary = task[(task.name.=="boundary") .& (task.epoch .== epoch), :]
append!(boundary, DataFrame(boundary[1,:]))

# Other data
#iso = Load.column_load_spikes("isolated", "RY16", 36)
#Load.register(iso, spikes, on="time", transfer=["isolated"])
Munge.spiking.isolated(spikes, lfp)

spikes[!,:isolated_area] = convert(Vector{Float16}, spikes.isolated)
spikes[spikes.isolated_area .== 1 .&& spikes.area.=="PFC",:isolated_area] .= 0.5


import Munge
Munge.behavior.annotate_relative_xtime!(beh)
cycles  = annotate_explodable_cycle_metrics(beh, cycles, dat, x, y, T)
ripples = annotate_explodable_cycle_metrics(beh, ripples, dat, x, y, T)
cells = Load.load_cells(animal, day, "*")
GC.gc()

# -------- END of non-checkpointed PREPROC ----------------------------

## FIGURE SETTINGS
visualize = :slider
# Spike colors
colorby    = "isolated" # nothing | a column of a dataframe
#colorby    = "vs(𝔾|ℙ,ℙ|𝔾)" # nothing | a column of a dataframe
#colorby    = "vs(ℙₛ|ℙ,ℙ|ℙₛ)" # nothing | a column of a dataframe
colorwhere = :spikes # dataframe sampling from
plotNon = false
#theta_phase_color = :specific_color_of_phase # how the theta WAVE is drawn
theta_phase_color = :specific_color_of_phase

# Color any sorting along yaxis
#colorby    = "bestTau_marginal=x-y_datacut=:all" # nothing | a column of a dataframe
#resort_cell_order = :meanrate
colorby    = "isolated" # nothing | a column of a dataframe
resort_cell_order = [:area, :meanrate]

# --------------------
# RUN VIDEO 
# --------------------
cells, spikes = Load.cell_resort(cells, spikes, resort_cell_order)

Δ_bounds = [0.20, 0.20] # seconds

start_index = Dict("beh"=>Utils.searchsortednearest(beh.time, T[1]),
          "lfp"=>Utils.searchsortednearest(lfp.time, T[1]))

# Period of sample
Δt = Dict("beh"=>median(diff(beh.time)),
          "lfp"=>median(diff(lfp.time)),
          "prob"=>median(diff(T)))

# Samples to step
Δi = Dict("beh"=>(Δt["beh"]/Δt["prob"])^-1,
          "lfp"=>(Δt["lfp"]/Δt["prob"])^-1)

function get_extrema(df, col) 
    x = dropmissing(spikes[!, [:isolated]])[!,colorby]
    x = Bool <: eltype(x) ? Int64.(x) : x
    nanextrema(x)
end


if resort_cell_order !== nothing
    cells, spikes = Load.cell_resort(cells, spikes, resort_cell_order)
end

## ---------------
## FIGURE WISHLIST
## ---------------
# - Goals
#   - Glow goal
#   - □ Vector to goal?
# - Errors
#   - ✔ Show actual poked well that caused the error
#   - □ Pass poke data over trial in seperate column
# - Sequences
#   - ✓ Vector: Start to end
#   - □ Vector: animal to end
#   - □ Scatter of maxima
# - Theta
#   □ Sweeps suck
#   □ θ sweep arrow : start at trough, end at peak
# - Ripples
# - Histograms for properties of interest
#   (I.e. best way to test for effects. Example, cosine similarity of animal vec to goal and replay vec)
# - Fields of cells
#   - Write package to generate convex hull regions of cells.
#   - Output to format for plots.jl and makie.jl
#   - (Option) show fields within theta cycle / ripple cycle
# - Events
#   - event :: number that counts across theta/ripple events
#   - (option) vspan current event or all events in plot window


# ======
# FIGURE
# ======
Fig = Figure(resolution=(800,1300))
if visualize == :video
    t = Observable(Int(1000))
elseif visualize == :slider
    sl_t = Slider(Fig[1:2, 2], range = 1:length(T), horizontal = false, 
                  startvalue = 1000)
    t = lift(sl_t.value) do tt
        tt
    end
end

# ======
# DATA
# ======
using Decode
using Decode.makie_observable
@time beh_window   = @lift select_range($t, T, data=beh, Δ_bounds=[-0.01, 0.60])
@time spike_window = @lift select_range($t, T, data=spikes, Δ_bounds=Δ_bounds)
@time lfp_window   = @lift select_est_range($t, T, start_index["lfp"], Δt["lfp"],
                                            data=lfp, Δ_bounds=Δ_bounds) # TODO warning this only works for an epoch alone, times skip are not constant across epoch borders
@time ripple_window = @lift Decode.select_events($t, T, events=ripples)
@time cycle_window  = @lift Decode.select_events($t, T, events=cycles)
theta_prob  = @lift select_prob($t,  T; prob=theta)
ripple_prob = @lift select_prob4($t, T; prob=ripple)

lfp_now    = @lift Int32(round(size($lfp_window,1)/2))
lfp_window = @lift transform($lfp_window, 
                             [:phase,:amp] => ((x,y)->x.*(y.^1.9)) => :phaseamp,
                             :time => (x->x.-mean(x)) => :time
                            )
now = 1
cmlabel = Dict(0=>"cue", 1=>"mem", -1=>"nontask")

# Grids
# Axes
TT = length(T)
correct = @lift $beh_window.correct[1]
# Arena data AXIS :: Behavior
cuemem_future  = @lift cmlabel[$beh_window.cuemem[now]]
@time title_str = @lift "i=$($t), %=$(@sprintf("%2.2f",100*$t/TT)) \ncorr=$($correct), goal=$($beh_window.stopWell[1]), vel=$(@sprintf("%2.1f",abs($beh_window.velVec[1])))"

axArena  = Axis(Fig[1,1], xlabel="x", ylabel="y", title=title_str)
gNeural = Fig[2,1] = GridLayout(1,1)
axLFPsum = Axis(gNeural[1,1], ylims=(50,150))
axLFP = Axis(gNeural[2:3,1])
axNeural = Axis(gNeural[3:8,1], xlabel="time")
controls = Fig[1:2, 3]

# TODO PROBLEMS

# THETA RAINBOW phase middle
if theta_phase_color == :rainbow_event
    start = @lift mean($lfp_window.time) .- Δ_bounds[1]
    stop  = @lift mean($lfp_window.time) .+ Δ_bounds[2]
    nan_color = RGBA(1,1,1,1)
    zerod_phase = @lift (($lfp_window.time .-$start))./($stop-$start)
    @lift $zerod_phase[($zerod_phase .> 1) .| ($zerod_phase .< 0)] .= NaN
    zerod_phase = get_color.(zerod_phase, :hawaii, nan_color)
    lines_theta = @lift lines!(axLFPsum, $lfp_window.time, $lfp_window.raw, colormap=:hawaii,
           color=$zerod_phase, linewidth=2)
elseif theta_phase_color == :rainbow
elseif theta_phase_color == :specific_color_of_phase
    ln_theta_color = @lift mean($lfp_window.color_phase[-10+$lfp_now:10+$lfp_now])
end
ln_phase_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_window[!,[:time,:phaseamp]])]) 
ln_theta    = lines!(axLFPsum,   ln_phase_xy, color=ln_theta_color, linestyle=:dash)
xlims!(axLFPsum, Δ_bounds .* (-1,1))

# -----------------
# LFP Raw data AXIS
# -----------------
ln_broad_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_window[!,[:time,:broadraw]])]) 
ln_theta_xy = @lift([Point2f(Tuple(x)) for x in eachrow($lfp_window[!,[:time,:raw]])]) 
ln_broad    = lines!(axLFP, ln_broad_xy, color=:gray,  linestyle=:dash)
ln_theta    = lines!(axLFP, ln_theta_xy, color=ln_theta_color, linestyle=:dash)
axLFP.yticks = 100:100
axLFPsum.yticks = 0:0
xlims!(axLFP, Δ_bounds .* (-1,1))

# ----------------
# Neural data AXIS
# ----------------
spike_colorrange = begin
    if colorby === nothing
        (-Inf, Inf)
    else
        if colorwhere == :behavior
            get_extrema(beh, colorby)
        elseif colorwhere == :spikes
            spikes[!, colorby] = disallowmissing(replace(spikes[!, colorby],
                                                         missing=>NaN))
            get_extrema(spikes, colorby)
        elseif colorwhere == :cells
            cells[!,colorby] = replace(cells[!, colorby], missing=>NaN)
            get_extrema(cells[spikes[!,:unit],:], colorby)
        else
            @error "Unrecognized symbol=$colorwhere"
        end
    end
end
spikes_abovebelow_quantile = nothing
if spikes_abovebelow_quantile !== nothing
    spike_colorrange=(-1, 1)
end
spike_colormap = if Utils.in_range([0], spike_colorrange)[1]
    top_value = max(abs.(spike_colorrange)...)
    spike_colorrange = (-top_value,top_value)
    (spike_cmap_with0,0.9)
else
    (spike_cmap_zeroless,0.9)
end

# COLORBAR :: TODO, make this derive from the spike_colors arg
#spike_cmap_with0 = :vik
#spike_cmap_zeroless = :viridis
#spike_cbar = colorby===nothing ? nothing : Colorbar(gNeural[3:8,2],
#                                                   limits=spike_colorrange,
#                                                   colormap=spike_cmap_zeroless,
#                                                   label=string(colorby),
#                                                   flipaxis=false,
#                                                   vertical=true)

spike_colors = @lift begin
    if colorby === nothing
        :white
    else
        if colorwhere == :behavior
            $beh_window[!, colorby]
        elseif colorwhere == :spikes
            $spike_window[!, colorby]
        elseif colorwhere == :cells
            cells[$spike_window[!, :unit], colorby]
        else
            @error "Unrecognized symbol=$colorwhere"
        end
    end
end

sc_sp_window = @lift([Point2f(Tuple(x)) for x in
                      eachrow($spike_window[!,[:time,:unit]])])

sc = scatter!(axNeural, sc_sp_window, colormap=spike_colormap, markersize=6,
              colorrange=spike_colorrange, color=spike_colors)
xlims!((Δ_bounds .* (-1,1))...)

if resort_cell_order == :meanrate
    hlines!(axNeural, findfirst(cells.meanrate .>6), color=:darkgrey, linestyle=:dash)
end

#ln_lfp_phase = @lift([Point2f(Tuple(x)) for x in
#                      eachrow(select($lfp_window, :time, :phase_plot=>x->x.*1))])
#sc = lines!(axNeural, ln_lfp_phase,  color=:gray, linestyle=:dash)
vlines!(axNeural, [0], color=:red, linestyle=:dash)
vlines!(axLFP,    [0], color=:red, linestyle=:dash)
vlines!(axLFPsum, [0], color=:red, linestyle=:dash)

# Arena data AXIS :: Probabilities
#hm_prob = @lift select_prob($t)
#hm = heatmap!(axArena, x,y, hm_prob, colormap=(:bamako,0.8))
ca1    = @lift($ripple_prob[:,:,1])
ca1pfc = @lift($ripple_prob[:,:,2])
pfcca1 = @lift($ripple_prob[:,:,3])
pfc    = @lift($ripple_prob[:,:,4])
limits = nanextrema(theta)
hm_theta     = heatmap!(axArena, x,y, theta_prob, colormap=(:romaO,0.8), colorrange=limits, interpolate=false)
hm_ripple    = heatmap!(axArena, x,y, ca1,        colormap=(:linear_ternary_blue_0_44_c57_n256, 0.8), interpolate=false)
hm_ripple_hp = heatmap!(axArena, x,y, ca1pfc,     colormap=(:summer, 0.8), interpolate=false)
hm_ripple_ph = heatmap!(axArena, x,y, pfcca1,     colormap=(:spring, 0.8), interpolate=false)
hm_ripple_p  = heatmap!(axArena, x,y, pfc,        colormap=(:linear_ternary_red_0_50_c52_n256, 0.8), interpolate=false)
if plotNon
    non_prob = @lift select_prob($t, T, prob=non)
    hm_non = heatmap!(axArena, x,y, non_prob, colormap=(:bamako,0.8))
end
xlims!(axArena, 40,175)
ylims!(axArena, 10, 95)


# ----------------
# Animal and wells
# ----------------
point_xy = @lift(Point2f($beh_window.x[now], $beh_window.y[now]))
line_xy  = @lift([Point2f(b.x, b.y) for b in eachrow($beh_window)])
sc_xy    = scatter!(axArena, point_xy, color=:white)
sc_arena = scatter!(axArena, arenaWells.x, arenaWells.y, marker=:star5, markersize=40, color=:gray)
sc_home  = scatter!(axArena, [homeWell.x], [homeWell.y],   marker='𝐇', markersize=25,   color=:gray)
ln_xy = lines!(axArena, line_xy, color=:white, linestyle=:dash)
lines!(axArena, boundary.x, boundary.y, color=:grey)

# -----------------
# Show place fields
# -----------------
if hasproperty(Main,:F)
end

# ----------------
# Sequence vectors
# ----------------
#cycle_start_xy  = @lift isempty($cycle_window) ? [NaN, NaN]  : $cycle_window.dec₀
#cycle_uv        = @lift isempty($cycle_window) ? [NaN, NaN]  : $cycle_window.dec₀₁
#cycle_arrows = @lift arrows!(axArena, real($cycle_start_xy[1]), imag($cycle_start_xy[1]),
#                             real($cycle_uv[1]), imag($cycle_uv[1]))
#
#ripple_start_xy  = @lift isempty($ripple_window) ? [NaN, NaN]  : $ripple_window.dec₀
#ripple_uv        = @lift isempty($ripple_window) ? [NaN, NaN]  : $ripple_window.dec₀₁
#ripple_arrows = @lift arrows!(axArena, real($ripple_start_xy[1]), imag($ripple_start_xy[1]),
#                              real($ripple_uv[1]), imag($ripple_uv[1]))

# ---------------
# Sequence maxima
# ---------------

# ----------------
# Goal vectors
# ----------------

# ---------------
# Goals and pokes
# ---------------
future₁ = @lift($beh_window.stopWell[now])
future₂ = @lift($beh_window.futureStopWell[now])
past₁   = @lift($beh_window.pastStopWell[now])
correctColor = @lift begin
    if $correct == 0
        :indianred1
    elseif $correct == 1
        :mediumspringgreen
    else
        :white
    end
end
future₁_well_xy = @lift $future₁ == -1 ? Point2f(NaN, NaN) : Point2f(wells.x[$future₁], wells.y[$future₁])
future₂_well_xy = @lift (($future₂ == -1) || ($future₂ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$future₂], wells.y[$future₂])
past₁_well_xy    = @lift (($past₁ == -1) || ($past₁ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$past₁], wells.y[$past₁])
sc_future₁_well = scatter!(axArena, future₁_well_xy,  alpha=0.5, marker='𝐅', markersize=60, color=correctColor, glowwidth=29)
sc_future₂_well = scatter!(axArena, future₂_well_xy,  alpha=0.5, marker='𝐟', markersize=60, color=:white,       glowwidth=5)
sc_past_well    = scatter!(axArena, past₁_well_xy,    alpha=0.5, marker='𝐏', markersize=60, color=:white,       glowwidth=5)
justify = 5
cuemem_future₂ = @lift cmlabel[$beh_window.futureCuemem[now]]
cuemem_past  = @lift cmlabel[$beh_window.pastCuemem[now]]
textsize=20
an_future₁_well = text!(axArena, cuemem_future, align=(:right,:bottom),  position=future₁_well_xy, 
                        textsize=textsize, color=:blue)
an_future₂_well = annotations!(axArena, cuemem_future₂,  future₂_well_xy, textsize, color=:blue )
an_past_well    = annotations!(axArena, cuemem_past,     past₁_well_xy, textsize, color=:blue  )
if doPrevPast
    past₂   = @lift($beh_window.pastStopWell[now])
    past₂_well_xy    = @lift (($past₂ == -1) || ($past₂ == $future₁)) ? Point2f(NaN, NaN) : Point2f(wells.x[$past₂], wells.y[$past₂])
end

# Pokes if error
poke_1_xy = @lift (!($correct==0) && $beh_window.poke_1[now]==1) ? Point2f(wells.x[1], wells.y[1]) : Point2f(NaN, NaN)
poke_2_xy = @lift (!($correct==0) && $beh_window.poke_2[now]==1) ? Point2f(wells.x[2], wells.y[2]) : Point2f(NaN, NaN)
poke_3_xy = @lift (!($correct==0) && $beh_window.poke_3[now]==1) ? Point2f(wells.x[3], wells.y[3]) : Point2f(NaN, NaN)
poke_4_xy = @lift (!($correct==0) && $beh_window.poke_4[now]==1) ? Point2f(wells.x[4], wells.y[4]) : Point2f(NaN, NaN)
poke_5_xy = @lift (!($correct==0) && $beh_window.poke_5[now]==1) ? Point2f(wells.x[5], wells.y[5]) : Point2f(NaN, NaN)
poke_1_sc = scatter!(axArena, poke_1_xy, marker='↑', markersize=30, color=:white, glow_width=10)
poke_2_sc = scatter!(axArena, poke_2_xy, marker='↑', markersize=30, color=:white, glow_width=10)
poke_3_sc = scatter!(axArena, poke_3_xy, marker='↓', markersize=30, color=:white, glow_width=10)
poke_4_sc = scatter!(axArena, poke_4_xy, marker='↑', markersize=30, color=:white, glow_width=10)
poke_5_sc = scatter!(axArena, poke_5_xy, marker='←', markersize=30, color=:white, glow_width=10)

function play_graphic(stop=length(T)-2000; framerate = 90)
    start = t[]
    timestamps = range(start, stop, step=1)
    P=Progress(stop-start, desc="Video")
    path = joinpath(dirname(decode_file), "decode_split=$(split_num)_start=$(start)_stop=$stop.mp4")
    record(Fig, path, timestamps; framerate=framerate, compression=10) do stamp
        try
            t[] = stamp
        catch
            @warn "timestamp failed with t=$(t[])"
        end
        next!(P)
    end
end

# -------------
# VISUALIZE
# -------------
if visualize == :video
    play_graphic()
elseif visualize == :slider
end

