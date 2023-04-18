include("prep2plot.jl")
DF1 = copy(DF1)

# ------------------------------------------------
# TRYING A SEPARATE APPROACH
# ------------------------------------------------

# ---------------
# 👉 CONSTANTS 👈
# ---------------
chunks = [:traj]
rows   = [:startWell, :stopWell]
yax    = [:startWell_tmpl, :stopWell_tmpl]
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
sort!(DF1, [:time, :areas])
match_cols = get_pmatch_cols([:startWell,:stopWell])
DF1[!,:pmatch] = all(Matrix(DF1[!, match_cols[1]]) .== 
                    Matrix(DF1[!, match_cols[2]]), dims=2) |> vec
GC.gc()
# Clean out any single time trajectories
DFc = groupby(DF1, [:traj])
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
muDFcc = combine(groupby(DFcc, [:areas, :i_tmpl, :time]), 
    [:value] => mean => :mean, [:value] => length => :n)
g1 = groupby(DFcc, [:areas, :i_tmpl])
g2 = groupby(muDFcc, [:areas, :i_tmpl])
K = intersect(keys(g1), keys(g2))
k = K |> first
for k in K
    k = NamedTuple(k)
    # Create mean cenetered values
    m, g = g1[k], g2[k]
    @assert m[:, :time] == g[:, :time]
    # Match on times
    g[!, :value] = g[!, :value] .- m[!, :mean]
end
DFcc = combine(g1, identity)
# ------------------------------------------------
# Combining templates?
# ------------------------------------------------
if template_combine != :none
    U = Dict(x=>unique(DFcc[!, x]) for x in propertynames(DFcc))
    g1 = groupby(DFcc, [:areas, :traj])
    for g in g1
        # Combine all non-active templates into a single template
        # called other
        if template_combine == :mean
            # Which template is active this :traj?    
            active = g[1, :i_test]
            # Which templates are not active?
            not_active = setdiff(U[:i_tmpl], active)
            # Combine all non-active templates into a single template
            # called other
            g[!, :i_tmpl] = g[!, :i_tmpl] .∈ (not_active,) ? 0 : g[!, :i_tmpl]
        elseif template_combine == :start_g_stop
            active_stop = g[1, :i_stop]
            # Set any points that are not the active stop to -1
            g[!, :i_tmpl] = g[!, :i_tmpl] .∉ (active_stop,) ? -1 : g[!, :i_tmpl]
        elseif template_combine == :stop_g_start
            active_start = g[1, :i_start]
            # Set any points that are not the active start to -1
            g[!, :i_tmpl] = g[!, :i_tmpl] .∉ (active_start,) ? -1 : g[!, :i_tmpl]
        end
    end
    DFcc = combine(g1, identity)
    # Remove any -1 :i_tmpl
    DFcc = @subset(DFcc, :i_tmpl .> 0)
    # Combine any time points that have the same :i_tmpl, :area, :i_test
    other_names = setdiff(names(DFcc), [:areas, :i_tmpl, :time, :value])
    DFcc = combine(groubpy(DFcc, [:areas, :i_tmpl, :time]), 
                            other_names .=> first .=> other_names,
                                  :value => length => :n,
                                  :value => mean => :value
    )
end
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡
# ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# ------------------------------------------------
# PLOT: Trajectories
# ------------------------------------------------
C = OrderedDict()
# Plotting columnar chunks
Chunks = groupby(means_subtracted ? DFcc : DFc, chunks)
chunk = Chunks |> first
# for (c,chunk) in enumerate(Chunks)
    Rows = groupby(chunk, rows)    
    R = OrderedDict()
    # Setup rows of a subplots
    (k, row) = zip(keys(Rows), Rows) |> first
    # for (k,row) in zip(keys(Rows), Rows)
        dur = extrema(row.time) |> x-> round(x[2] - x[1], digits=2)
        p = plot(title="dur=$dur")
        m = maximum(row.value)
        Yax = groupby(row, yax)
        i,y = 1,(Yax |> first)
        global startmove = stopmove = nothing
        for (i, y) in enumerate(Yax)
            plot!(p, y.time, (i-1).*m .+ y.value.*mscale; 
                   markersize=1,
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
            actual == label ? "ACTUAL: $label" : label
        end
        global startmove, stopmove
        # vspan!(p, startmove, stopmove; fillalpha=0.1, linealpha=0.2) # BUG: this does not match the record below
        ytcks = ((axes(Yax, 1) .- 1) .* m, ylabels)
        actual = "$(row.startWell[1])-$(row.stopWell[1])"
        plot!(;yticks=ytcks, xlabel="time", ylabel="react per tmpl", 
            title="dur=$dur, act=$actual", legend=false)
        traj = row.traj[1]
        trajplot=subset(beh, :traj => t->t.== traj, view=true) |> 
            @df plot(:time, [:speed], ylim=(0,40),legend=false,
                fill=0, fillstyle=:\, fillalpha=0.2) 
        plot!(y.time, y.moving.*15, ylim=(0,40),legend=false,
            fill=0, fillalpha=0.8, c=:red)
        layout = @layout [a0.9h; b{0.1h}]
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
