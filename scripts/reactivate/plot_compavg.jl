include("prep2plot.jl")
DF1 = copy(DF1)

# ------------------------------------------------
# TRYING A SEPARATE APPROACH
# ------------------------------------------------

# ---------------
# 👉 CONSTANTS 👈
# ---------------
g_chunk = [:traj]
g_rows   = [:startWell, :stopWell]
g_yax    = [:startWell_tmpl, :stopWell_tmpl]
subs =   [:areas => a-> a .== "ca1-ca1", 
          :moving_tmpl => a-> a .== true]
mscale = 8
means_subtracted = true
template_combine = :mean # :none, :mean, :start_g_stop, stop_g_start
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
        [:value] => mean => :mean, [:value] => length => :n)
    g1 = groupby(DFcc, [:areas, :i_tmpl])
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
            g[!, :value] = g[!, :value] .- m[!, :mean]
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
    println("Rows before combine process: $(nrow(DFcc))")
    # BUG: i_test for every traj has 2 values, is not unique!
    @info "Combining templates using $template_combine"
    U = DFcc |> uniq
    g1 = groupby(DFcc, [:areas, :traj])
    g = g1[1]
    for g in g1
        # Combine all non-active templates into a single template
        # called other
        if template_combine == :mean
            # Which template is active this :traj?    
            active = unique(g[:, :i_test])
            # Which templates are not active?
            not_active = setdiff(U[:i_tmpl], active)
        # Combine all non-active templates into a single template, except
        # for the reverse template. All others lumped into "other"
        elseif template_combo  ==  :mean_leavereverse
            # Which template is active this :traj?    
            active = unique(g[:, :i_test])
            startWell, stopWell = unique(g[:, :startWell]), unique(g[:, :stopWell])
            i_tmpl_reverse = subset(g, :startWell_tmpl => x -> x == stopWell, 
                :stopWell_tmpl => x -> x == startWell, view=true).i_tmpl |> 
                unique 
            # Which templates are not active?
            not_active = setdiff(U[:i_tmpl], active)
        elseif template_combine == :start_g_stop
            active_stop = g[1, :i_stop]
            # Set any points that are not the active stop to -1
            g[!, :i_tmpl] = g[!, :i_tmpl] .∉ (active_stop,) ? -1 : g[!, :i_tmpl]
        elseif template_combine == :stop_g_start
            active_start = g[1, :i_start]
            # Set any points that are not the active start to -1
            g[!, :i_tmpl] = g[!, :i_tmpl] .∉ (active_start,) ? -1 : g[!, :i_tmpl]
        end
        if template_combine in [:mean, :means_subtracted]
            # Combine all non-active templates into a single template
            # called other
            f(x) = x ∈ not_active ? 0 : x
            g[!, :i_tmpl] = f.(g[!, :i_tmpl])
            g[!, :startWell_tmpl][g[!, :i_tmpl] .== 0] .= 0
            g[!, :stopWell_tmpl][g[!, :i_tmpl] .== 0]  .= 0
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
    println("Rows after combine process: $(nrow(DFcc))")
end
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# ------------------------------------------------
# PLOT: Trajectories
# ------------------------------------------------
C = OrderedDict()
# Plotting columnar chunks
Chunks = groupby(means_subtracted ? DFcc : DFc, g_chunk)
(c, chunk) = (1, Chunks |> first)
# for (c,chunk) in enumerate(Chunks)
    Rows = groupby(chunk, g_rows)    
    R = OrderedDict()
    # Setup rows of a subplots
    (k, row) = zip(keys(Rows), Rows) |> first
    # for (k,row) in zip(keys(Rows), Rows)
        dur = extrema(row.time) |> x-> round(x[2] - x[1], digits=2)
        p = plot(title="dur=$dur")
        m = maximum(row.value) # BUG: 0 when means_subtracted
        Yax = groupby(row, g_yax)
        i,y = 1,(Yax |> first)
        global startmove = stopmove = nothing
        for (i, y) in enumerate(Yax)
            # @infiltrate
            plot!(p, y.time, (i-1).*m .+ y.value.*mscale; 
                   markersize=1, fill=
                fillalpha=0.5, linealpha=0.2, legend=false)
            yc = y
            difs  = diff([Int8(0); Int8.(yc.moving)])
            global startmove = yc[findall(difs .== 1), :time]
            global stopmove  = yc[findall(difs .== -1), :time]
        end
        ylabels = map(Yax |> collect) do y
            y=first(eachrow(y))
            actual = "$(y.startWell)-$(y.stopWell)"
            label  = "$(y.startWell_tmpl)-$(y.stopWell_tmpl)"
            if y.i_tmpl == 0 
                "OTHER"
            else
                actual == label ? "ACTUAL: $label" : label
            end
        end
        global startmove, stopmove
        # vspan!(p, startmove, stopmove; fillalpha=0.1, linealpha=0.2) # BUG: this does not match the record below
        ytcks = (((axes(Yax, 1)|>collect) .- 1) .* m, ylabels)
        actual = "$(row.startWell[1])-$(row.stopWell[1])"
        plot!(;yticks=ytcks, xlabel="time", ylabel="react per tmpl", 
            title="dur=$dur, act=$actual", legend=false)
        traj = row.traj[1]
        trajplot=subset(beh, :traj => t->t.== traj, view=true) |> 
            @df plot(:time, [:speed], ylim=(0,40),legend=false,
                fill=0, fillstyle=:\, fillalpha=0.2) 
        plot!(y.time, y.moving.*15, ylim=(0,40),legend=false,
            fill=0, fillalpha=0.8, c=:red)
        layout = @layout [a{0.9h}; b{0.1h}]
        R[k] = plot(p,trajplot,layout=layout)
    # end
    C[c] = R
# end


# ------------------------------------------------
# PLOT: Corrplot of the variables above
# ------------------------------------------------
# DFr = DFc[rand(1:nrow(DFc),10_000), :]
# @df DFr corrplot([:startWell :stopWell :startWell_tmpl :stopWell_tmpl], 
#     grid = false, method = :pearson, order = :hclust)
