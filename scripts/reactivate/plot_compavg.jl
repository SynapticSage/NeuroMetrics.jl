include("prep2plot.jl")

# ------------------------------------------------
# TRYING A SEPARATE APPROACH
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

# SETTINGS
subs =   [:areas => a-> a .== "ca1-ca1", 
          :moving_tmpl => a-> a .== true,
    :traj => t -> (T=(t .∉ (getindex.(remove,1),)) )
]
DFc = subset(DF1, subs..., view=true)
DFc.value = DFc.mean
DFcc = copy(DFc)
sort!(DFc, [:areas, :time]) 

# CONSTANTS
chunks = [:traj]
rows   = [:startWell, :stopWell]
yax    = [:startWell_tmpl, :stopWell_tmpl]
mscale = 8

# INFO: TEST OF LOOP
C = OrderedDict()
# Plotting columnar chunks
Chunks = groupby(DFc, chunks)
chunk = Chunks |> first
for (c,chunk) in enumerate(Chunks)
    Rows = groupby(chunk, rows)    
    R = OrderedDict()
    # Setup rows of a subplots
    (k, row) = zip(keys(Rows), Rows) |> first
    for (k,row) in zip(keys(Rows), Rows)
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
        vspan!(p, startmove, stopmove; fillalpha=0.1, linealpha=0.2)
        yticks = ((axes(Yax, 1) .- 1) .* m, ylabels)
        actual = "$(row.startWell[1])-$(row.stopWell[1])"
        plot!(;yticks, xlabel="time", ylabel="react per tmpl", 
            title="dur=$dur, act=$actual", legend=false)
        traj = row.traj[1]
        trajplot=subset(beh, :traj => t->t.== traj, view=true) |> 
            @df plot(:time, [:speed :moving.*15], ylim=(0,25),legend=false) 
        plot!(y.time, y.moving.*15, ylim=(0,25),legend=false)
        layout = @layout [a{0.9h}; b{0.1h}]
        R[k] = plot(p,trajplot,layout=layout)
    end
    C[c] = R
end


# ------------------------------------------------
# PLOT: Corrplot of the variables above
# ------------------------------------------------
# DFr = DFc[rand(1:nrow(DFc),10_000), :]
# @df DFr corrplot([:startWell :stopWell :startWell_tmpl :stopWell_tmpl], 
#     grid = false, method = :pearson, order = :hclust)
