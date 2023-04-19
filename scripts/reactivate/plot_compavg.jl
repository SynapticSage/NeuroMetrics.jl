include("prep2plot.jl")
using GoalFetchAnalysis.Plot
Plot.setparentfolder("reactivation", "compavg")
Plot.setappend(opt)

# Just so language-server protocol can find the symbols
# from the script that generated the data.
if !isdefined(Main, :DF1); load_react_vars(["DF1"]); end
if !isdefined(Main, :beh); beh=DI.load_behavior(opt["animal"], opt["day"]); end
DIutils.pushover("Finished loading reactivate variables")

# |------------------------------------------------|
# |-----------  TRIAL-WISE PLOTTING ---------------|
# |------------------------------------------------|

# ||---------------||
# ||👉 CONSTANTS 👈||
# ||---------------||
# Plot control arguments
# ----------------------
# - `g_chunk` : chunks of plots to make
# - `g_rows`  : rows of plots to make
# - `g_yax`   : split vertically on y-axis of plots to make
# - `split_inner_loop` : split plotting within the y-axis loop
# - `subs`    : subset the data before plotting
# - `sumfunc` : function to use for summarizing data
# - `mscale`  : scale factor for mean centering vertical groups on y-axis
# - `means_subtracted` : whether to subtract the mean of each group
# - `template_combine` : how to combine templates
g_chunk  = [:traj]
g_rows   = [:startWell, :stopWell]
g_yax    = [:startWell_tmpl, :stopWell_tmpl]
split_inner_loop = [:traj, :startWell_tmpl, :stopWell_tmpl, :moving]
subs =   [:areas => a-> a .== "ca1-ca1", 
          :moving_tmpl => a-> a .== true]
sumfunc = mean # mean | median: used for mean centering
mscale  = 8
means_subtracted = true
template_combine = :mean # :none, :mean, :start_g_stop, stop_g_start
# ------------------------------------------------
# template_combine descriptions
# ------------------------------------------------
# - none : no combining
# - mean : mean all non-active templates(ie tmpl=test) 
#          into a single template called other
# - start_g_stop: keep only the start tmpl given the actual stop template
# - stop_g_start: keep only the stop tmpl given the actual start template
# ------------------------------------------------
sort(DF1, [:time, :areas])
match_cols = get_pmatch_cols([:startWell,:stopWell])
DF1[!,:pmatch] = all(Matrix(DF1[!, match_cols[1]]) .== 
                     Matrix(DF1[!, match_cols[2]]), dims=2) |> vec
GC.gc()
# Clean out any single time trajectories
DFc = groupby(copy(DF1), [:traj])
remove = []
for (k, df) in zip(keys(DFc), DFc)
    if length(unique(df.time)) ≤ 1
        push!(remove, k)
    end
end
println("Removing $(length(remove)) single time trajectories")
push!(subs, :traj => t -> (T=(t .∉ (getindex.(remove,1),)) ))
DFc = subset(DF1, subs..., view=true)
DFc.value = DFc.mean
sort!(DFc, [:areas, :time]) 

# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡ APPLY TEMPATE COMBINE AND MEAN CENTERING ≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# Create a mean centered version of the data
DFcc = copy(DFc)
if means_subtracted
    @info "Mean centering data"
    muDFcc = combine(groupby(DFcc, [:areas, :time]), 
        [:value] => sumfunc => :value, [:value] => length => :n)
    g1 = groupby(DFcc,   [:areas, :i_tmpl])
    g2 = groupby(muDFcc, [:areas])
    K = keys(g1)
    k = K |> first
    for k in K
        k = NamedTuple(k)
        k2 = (;areas=k.areas,)
        # Create mean cenetered values
        g, m = g1[k], g2[k2]
        @assert m[:, :time] == g[:, :time]
        # Match on times
        try
            g[!, :value] = g[!, :value] .- m[!, :value]
        catch
            @infiltrate
        end
    end
    DFcc = combine(g1, identity)
end

# ------------------------------------------------
# Combining templates?
# ------------------------------------------------
uniq(DFcc) = Dict(x=>unique(DFcc[!, x]) for x in propertynames(DFcc))
if template_combine != :none
    rowsbefore = nrow(DFcc);
    @info "Combining templates using $template_combine"
    U = DFcc |> uniq
    g1 = groupby(DFcc, [:areas, :traj])
    g = g1[1]
    @showprogress for g in g1
        # Combine all non-active templates into a single template
        # called other
        try
            if template_combine == :mean
                # Which template is active this :traj?    
                active = unique(g[:, :i_test])
                # Which templates are not active?
                not_active = setdiff(U[:i_tmpl], active)
            # Combine all non-active templates into a single template, except
            # for the reverse template. All others lumped into "other"
            elseif template_combine  ==  :mean_w_reverse
                # Which template is active this :traj?    
                active = unique(g[!, :i_test])
                startWell, stopWell = unique(g[!, :startWell]), unique(g[!, :stopWell])
                if length(startWell) > 1 || length(stopWell) > 1
                    @warn "traj $g[1,:traj]: More than one startWell or stopWell found"
                    i = maximum(countmap(g[!, :startWell])) |> first
                    j = maximum(countmap(g[!, :stopWell])) |> first
                    startWell, stopWell = i, j
                end
                i_tmpl_reverse = subset(g, :startWell_tmpl => x -> x .== stopWell, 
                    :stopWell_tmpl => x -> x .== startWell, view=true).i_tmpl |> 
                    unique 
                # Which templates are not active?
                not_active = setdiff(U[:i_tmpl], [active..., i_tmpl_reverse...])
            elseif template_combine == :start_g_stop
                active_stop = g[1, :i_stop]
                # Set any points that are not the active stop to -1
                g[!, :i_tmpl] = g[!, :i_tmpl] .∉ (active_stop,) ? -1 : g[!, :i_tmpl]
            elseif template_combine == :stop_g_start
                active_start = g[1, :i_start]
                # Set any points that are not the active start to -1
                g[!, :i_tmpl] = g[!, :i_tmpl] .∉ (active_start,) ? -1 : g[!, :i_tmpl]
            end
            if template_combine in [:mean, :mean_w_reverse]
                # Combine all non-active templates into a single template
                # called other
                f(x) = x ∈ not_active ? 0 : x
                new_tmpl = g[!, :i_tmpl] = f.(g[!, :i_tmpl])
                g[!, :startWell_tmpl][new_tmpl .== 0] .= 0
                g[!, :stopWell_tmpl][new_tmpl .== 0]  .= 0
            end
        catch e
            @warn "Error $e in traj $(g[1,:traj])"
        end
    end
    DFcc = combine(g1, identity)
    # Remove any -1 :i_tmpl
    DFcc = @subset(DFcc, :i_tmpl .>= 0)
    # Combine any time points that have the same :i_tmpl, :area, :i_test
    other_names = setdiff(names(DFcc), ["areas", "i_tmpl", "time", "value"])
    DFcc = combine(
        groupby(DFcc, [:areas, :i_tmpl, :time]), 
                            other_names .=> first .=> other_names,
                            :value => length => :nn,
                            :value => mean => :value
    )
    println("Rows before combine process: $(rowsbefore)")
    println("Rows after  combine process: $(nrow(DFcc))")
end
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# ------------------------------------------------
# PLOT: - - - - - Trajectories  - - - - - - - - -
# ------------------------------------------------
Plot.setfolder("trajectories")
C = OrderedDict()
# Plotting columnar chunks
Chunks = groupby(means_subtracted ? DFcc : DFc, g_chunk)
(c, chunk) = (1, Chunks |> first)
for (c,chunk) in enumerate(Chunks)
    Rows = groupby(chunk, g_rows)    
    R = OrderedDict()
    # Setup rows of a subplots
    (k, row) = zip(keys(Rows), Rows) |> first
    for (k,row) in zip(keys(Rows), Rows)
        dur = extrema(row.time) |> x-> round(x[2] - x[1], digits=2)
        p = plot(title="dur=$dur")
        m = maximum(row.value) 
        Yax = groupby(row, g_yax)
        i,y = 1,(Yax |> first)
        global startmove = stopmove = nothing
        for (i, y) in enumerate(Yax)
            # @infiltrate
            plot!(p, y.time, (i-1).*m .+ y.value.*mscale; 
                markersize=1, fill=((i-1).*m),
                fillalpha=0.5, linealpha=0.2, legend=false)
            yc = y
            difs  = diff([Int8(0); Int8.(yc.moving)])
            global startmove = yc[findall(difs .== 1), :time]
            global stopmove  = yc[findall(difs .== -1), :time]
        end
        global startmove, stopmove
        # vspan!(p, startmove, stopmove; fillalpha=0.1, linealpha=0.2) # BUG: this does not match the record below
        ylabels = tmpl_ylabels(Yax)
        ytcks = (((axes(Yax, 1)|>collect) .- 1) .* m, ylabels)
        actual = "$(row.startWell[1])-$(row.stopWell[1])"
        plot!(;yticks=ytcks, xlabel="time", ylabel="react per tmpl", 
            title="dur=$dur, act=$actual", legend=false)
        traj = row.traj[1]
        trajplot=subset(beh, :traj => t->t.== traj, view=true) |> 
            @df plot(:time, [:speed], ylim=(0,40),legend=false,
                fill=0, fillstyle=corerr_fillstyle(row.correct[1]), 
                fillalpha=0.2) 
        plot!(y.time, y.moving.*15, ylim=(0,40),legend=false,
            fill=0, fillalpha=0.8, c=:red)
        layout = @layout [a{0.9h}; b{0.1h}]
        R[k] = plot(p,trajplot,layout=layout)
    end
    # if multiple rows of plot, store the dict of rows, else just store 
    # the plot
    C[c] = length(R) > 1 ? R : R |> values |> first
end

# ------------------------------------------------
# Do the above with CA1-CA1 and CA1-PFC
# ------------------------------------------------


# ----------------------------------------------------------------
# Count (+ : above mean) and (- : below mean) reactivation events 
# during immobility before and during movement
# ----------------------------------------------------------------
Plot.setfolder("summary")
import Peaks
if means_subtracted
    g = groupby(means_subtracted ? DFcc : DFc, 
        unique([:i_tmpl, g_chunk..., g_rows..., g_yax..., :moving]))
    """
    Returns
    -------
    - p: peak prominences
    - w: peak widths
    - m: peak maxima
    """
    function get_peaks(x)
        if length(x) < 3
            @warn "Not enough data points to find peaks"
            return (; p=[], w=[], m=[], v=[])
        end
        println("Finding peaks")
        pks, m = Peaks.findmaxima(x)
        if isempty(pks)
            @warn "No peaks found"
            return (; p=[], w=[], m=[], v=[])
        end
        _, p = Peaks.peakproms(pks, x)
        _, w, _, _ = Peaks.peakwidths(pks, x, p)
        v = m .* w
        return (; p, w, m, v)
    end
    m   = combine(g,  :value => get_peaks => [:p, :w, :m, :v])
    mi  = combine(g, :value => (x->get_peaks(-x)) => [:p, :w, :m, :v])
    transform!(m, :i_tmpl => (x->x .> 0) => :i_tmpl)
    transform!(mi, :i_tmpl => (x->x .> 0) => :i_tmpl)
    gm  = groupby(m,  [:traj,  :i_tmpl,:moving])
    gmi = groupby(mi, [:traj, :i_tmpl, :moving])
    global gm  = combine(gm, :p => mean => :p_mean, :p => std => :p_std, :w =>
        mean => :w_mean, :w => std => :w_std, :m => mean => :m_mean, :m => std
            => :m_std, :p => length => :nrow, :v => mean => :v_mean, :v => std
                => :v_std)
    global gmi  = combine(gmi, :p => mean => :p_mean, :p => std => :p_std, :w
        => mean => :w_mean, :w => std => :w_std, :m => mean => :m_mean, :m =>
            std => :m_std, :p => length => :nrow, :v => mean => :v_mean, :v =>
                std => :v_std)

    # PLOT: i_tmpl vs p_mean +/- p_std, peak prominence
    using Measurements
    gmc = combine(groupby(subset(gm,:moving=>x->x.==false) , [:i_tmpl, :moving]),
                    :p_mean => length => :count, :p_mean => mean => :p_mean,
    [:p_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :p_std
    )
    gmci = combine(groupby(subset(gmi, :moving=>x->x.==false), [:i_tmpl, :moving]),
                    :p_mean => length => :count, :p_mean => mean => :p_mean,
    [:p_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :p_std
    )
    b1=@df gmc bar(:i_tmpl, :p_mean .± :p_std; xlabel="i_tmpl", 
        ylabel="Peak Heights", title="Mean  of (+)\nreactivation events", 
        legend=:bottomright, size=(300, 600), left_margin=10Plots.mm, 
        label=false, marker=:none, markeralpha=0.0, msw=5, 
        ylims=(0,0.1))
    b2= @df gmci bar(:i_tmpl, :p_mean .± :p_std, xlabel="i_tmpl", ylabel="Trough Heights",
        title="Mean  of (-)\nreactivation events", legend=:bottomright,
        size=(300, 600), left_margin=10Plots.mm, label=false, msw=5,
        yflip=true, ylims=(-0, -0.1), c=:red)
    ylims!(0,0.1)
    plot(b1,b2, layout=(1,2), size=(450, 600))
    Plot.save("PROMINENCE, peak-trough, stillness")

    # PLOT: i_tmpl vs m_mean, maximum of peak
    gmc = combine(groupby(subset(gm,:moving=>x->x.==false) , [:i_tmpl, :moving]),
                    :m_mean => length => :count, :m_mean => mean => :m_mean,
    [:m_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :m_std
    )
    gmci = combine(groupby(subset(gmi, :moving=>x->x.==false), [:i_tmpl, :moving]),
                    :m_mean => length => :count, :m_mean => mean => :m_mean,
    [:m_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :m_std
    )
    b1=@df gmc bar(:i_tmpl, :m_mean .± :m_std; xlabel="i_tmpl", 
        ylabel="Peak Heights", title="Mean  of (+)\nreactivation events", 
        legend=:bottomright, size=(300, 600), left_margin=10Plots.mm, 
        label=false, marker=:none, markeralpha=0.0, msw=5, 
        ylims=(0,0.1))
    b2= @df gmci bar(:i_tmpl, :m_mean .± :m_std, xlabel="i_tmpl", ylabel="Trough Heights",
        title="Mean  of (-)\nreactivation events", legend=:bottomright,
        size=(300, 600), left_margin=10Plots.mm, label=false, msw=5,
        yflip=true, ylims=(-0, -0.1), c=:red)
    ylims!(0,0.1)
    plot(b1,b2, layout=(1,2), size=(450, 600))
    Plot.save("MAXIMUM, peak-trough, stillness")

    # PLOT: i_tmpl vs w_mean, widths of peaks
    gmc = combine(groupby(subset(gm,:moving=>x->x.==false) , [:i_tmpl, :moving]),
                    :w_mean => length => :count, :w_mean => mean => :w_mean,
    [:w_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :w_std
    )
    gmci = combine(groupby(subset(gmi, :moving=>x->x.==false), [:i_tmpl, :moving]),
                    :w_mean => length => :count, :w_mean => mean => :w_mean,
    [:w_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :w_std
    )
    b1=@df gmc bar(:i_tmpl, :w_mean .± :w_std; xlabel="i_tmpl", 
        ylabel="Peak Heights", title="Mean  of (+)\nreactivation events", 
        legend=:bottomright, size=(300, 600), left_margin=10Plots.mm, 
        label=false, marker=:none, markeralpha=0.0, msw=5, 
        ylims=(0,0.1))
    b2= @df gmci bar(:i_tmpl, :w_mean .± :w_std, xlabel="i_tmpl", ylabel="Trough Heights",
        title="Mean  of (-)\nreactivation events", legend=:bottomright,
        size=(300, 600), left_margin=10Plots.mm, label=false, msw=5,
        yflip=true, ylims=(-0, -0.1), c=:red)
    Plot.save("WIDTH, peak-trough, stillness")

    # PLOT: Volumemetric means
    gmc = combine(groupby(subset(gm,:moving=>x->x.==false) , [:i_tmpl, :moving]),
                    :v_mean => length => :count, :v_mean => mean => :v_mean,
    [:v_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :v_std
    )
    gmci = combine(groupby(subset(gmi, :moving=>x->x.==false), [:i_tmpl, :moving]),
                    :v_mean => length => :count, :v_mean => mean => :v_mean,
    [:v_std,:nrow] => ((x,y)->sqrt(nansum(x.^2))/length(y)) => :v_std
    )
    b1=@df gmc bar(:i_tmpl, :v_mean .± :v_std; xlabel="i_tmpl", 
        ylabel="Peak Heights", title="Mean  of (+)\nreactivation events", 
        legend=:bottomright, size=(300, 600), left_margin=10Plots.mm, 
        label=false, marker=:none, markeralpha=0.0, msw=5, 
        )
    b2= @df gmci bar(:i_tmpl, :v_mean .± :v_std, xlabel="i_tmpl", ylabel="Trough Heights",
        title="Mean  of (-)\nreactivation events", legend=:bottomright,
        size=(300, 600), left_margin=10Plots.mm, label=false, msw=5,
        yflip=true, c=:red)
    plot(b1,b2, layout=(1,2), size=(450, 600))
    Plot.save("VOLUME, peak-trough, stillness")

end

# ------------------------------------------------
# PLOT: Characteristic frequency scale,
# FFT of each i_tmpl at each time
# ------------------------------------------------
Plot.setfolder("fft")
using FFTW, DSP
HH = []
P = []
for move in [true, false]
    tmp = subset(DFcc, :moving=>x->x.==move)
    movestr = move ? "moving" : "still"
    H = []
    for i1 in groupby(tmp, :i_tmpl)
        issorted(i1.time) || error("time not sorted")
        fs=1/median(diff(i1.time))
        S = DSP.mt_spectrogram(i1.value, 30; fs)
        h=plot(S.time, S.freq, S.power, xlabel="Time(s)", ylabel="Freq(Hz)",
        title="")
        ylims!(0,7.5)
        # ff = rfft(DFcc.value)
        # N = length(ff)
        # FFTW.
        push!(H, h)
    end
    b=Plot.blank(title="\n\n\n\nSpectrogram of of $movestr")
    p= plot(H[1:10]..., b, size=(1100, 500))
    Plot.save("fft-$movestr.png")
    push!(HH, H)
    push!(P, p)
end

# ------------------------------------------------
# PLOT: Corrplot of the template reponses
# ------------------------------------------------
# Show each test best active for its own i_tmpl
M1    = unstack(DFcc, :i_test, :i_tmpl, :value, combine=nanmean)
xtcks = tmpl_labels_dict(DFcc, "tmpl")
this_xtcks = OrderedDict(i=>xtcks[i] for i in tryparse.(Int, names(M1))
        if i in keys(xtcks))
_xtcks = (1:length(this_xtcks)|>collect, values(this_xtcks) |> collect)
ytcks = tmpl_labels_dict(DFcc, "test")
this_ytcks = OrderedDict(i=>ytcks[i] for i in M1[!, :i_test])
_ytcks = (1:length(this_ytcks)|>collect, values(this_ytcks) |> collect)
heatmap(Matrix((M1)[:, Not(:i_test)]), xlabel="i_tmpl", ylabel="i_test", 
    title="Mean of reactivation events", legend=:none, size=(600, 600),
    xticks=_xtcks, yticks=_ytcks, xrotation=90, yrotation=0)

M2=sort((unstack(DFc, :i_test, :i_tmpl, :value, combine=mean)), :i_test)
xtcks = tmpl_labels_dict(DFcc, "tmpl")
this_xtcks = OrderedDict(i=>xtcks[i] for i in tryparse.(Int, names(M2))
        if i in keys(xtcks))
_xtcks = (1:length(this_xtcks)|>collect, values(this_xtcks) |> collect)
ytcks = tmpl_labels_dict(DFcc, "test")
this_ytcks = OrderedDict(i=>ytcks[i] for i in M2[!, :i_test])
_ytcks = (1:length(this_ytcks)|>collect, values(this_ytcks) |> collect)
# Q: Those with the least off-diagonal correlation = better performance?
xi, yi = sortperm(_xtcks[2]), sortperm(_ytcks[2])
_xtcks = (_xtcks[1], _xtcks[2][xi])
_ytcks = (_ytcks[1], _ytcks[2][yi])
M2 = M2[:,Not("i_test")][yi, xi]
heatmap(Matrix(M2), xlabel="i_tmpl", ylabel="i_test", 
    title="Mean of reactivation events", legend=:none, size=(600, 600),
    xticks=_xtcks, yticks=_ytcks, xrotation=90, yrotation=0)


# ------------------------------------------------
# DFr = DFc[rand(1:nrow(DFc),10_000), :]
# @df DFr corrplot([:startWell :stopWell :startWell_tmpl :stopWell_tmpl], 
#     grid = false, method = :pearson, order = :hclust)
#
