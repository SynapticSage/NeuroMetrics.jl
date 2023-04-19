include("prep2plot.jl")
# ------------------------------------------------
# TRYING A SEPARATE APPROACH
# ------------------------------------------------
sort!(DF, [:time, :areas, :component])
match_cols = get_pmatch_cols([:startWell,:stopWell])
DF[!,:pmatch] = all(Matrix(DF[!, match_cols[1]]) .== 
                    Matrix(DF[!, match_cols[2]]), dims=2) |> vec
GC.gc()

# Clean out any single time trajectories
DFc = groupby(DF, [:traj])
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
          :traj => t->t ∉ Tuple(getindex.(remove,1))
]
chunks = [:traj]
rows   = [:startWell, :stopWell]
yax    = [:startWell_tmpl, :stopWell_tmpl]
DFc = subset(DF, subs..., view=true)

# BUG: losing compentents here!
# Step 2
push!(subs, :component => a-> a .<= 5)
DFc = subset(DF, subs..., view=true)
sort!(DFc, [:areas, :component, :time]) # BUG: α how does sorting affect this?
DIutils.pushover("Ready to go")

# Question: How many template combos does each trajectory have (it should be
# stable or nearly all)?
@time tmp=combine(groupby(DF, [:areas, :traj, :ha, :moving_tmpl]),
[:startWell_tmpl, :stopWell_tmpl] => 
    ((x,y) -> length(unique(eachrow([x y])))) => :tmpl_combos)
h1=histogram(tmp.tmpl_combos, group=tmp.areas, bins=1:1:20, normed=true, 
    bar_position=:stack,
    alpha=0.2,
    title="Histogram of template combos per trajectory", 
    xlabel="Number of template combos", ylabel="Frequency")

# And for the DFc?
@time tmp=combine(groupby(DFc, [:areas, :traj, :ha, :moving_tmpl]),
[:startWell_tmpl, :stopWell_tmpl] => 
    ((x,y) -> length(unique(eachrow([x y])))) => :tmpl_combos)
h1=histogram(tmp.tmpl_combos, group=tmp.areas, bins=1:1:20, normed=true, 
    bar_position=:stack,
    alpha=0.2,
    title="Histogram of template combos per trajectory", 
    xlabel="Number of template combos", ylabel="Frequency")


sort!(DFc, [:areas, :component, :time]) # BUG: α how does sorting affect this?

# Create a mean centered version of the data
DFcc = copy(DFc)
muDFcc = combine(groupby(DFcc, [:areas, :i_tmpl]), 
    [:value] => mean => :mean, [:value] => std => :std)
for g in groupby(DFcc, [:areas, :i_tmpl])
    # Create mean cenetered values
    @assert muDFcc[g, :time] == DFcc[g, :time]
    DFcc[g, :value] = DFcc[g, :value] .- muDFcc[g, :mean]
end

# --------------------------------
# PLOT: TEST examine traj
# --------------------------------
mean_centered = true
C = OrderedDict()
group = :component
# Plotting columnar chunks
Chunks = groupby(mean_centered ? DFcc : DFc, chunks)
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
            plot!(p, y.time, (i-1).*m .+ y.value; 
                   group=y[!,group], markersize=1,
                fillalpha=0.2, linealpha=0.2, legend=false)
            yc = groupby(y, :component) |> first
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
        vspan!(p, startmove, stopmove; fillalpha=0.2, linealpha=0.2)
        yticks = ((axes(Yax, 1) .- 1) .* m, ylabels)
        actual = "$(row.startWell[1])-$(row.stopWell[1])"
        plot!(;yticks, xlabel="time", ylabel="react per tmpl", 
            title="dur=$dur, act=$actual", legend=false)
        # R[k] = p
    # end
    # C[c] = R
# end

# PLOT: Examine trajectories
C = OrderedDict()
group = :component
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
        y = Yax |> first
        for (i, y) in enumerate(Yax)
            plot!(p, y.time, (i-1)*m .+ y.value; 
                group=y[!,group], markersize=1,
                 # fill=0, 
                fillalpha=0.2, linealpha=0.2, legend=false)
            yc = groupby(y, :component) |> first
            difs  = diff([Int8(0); Int8.(yc.moving)])
            startmove = findall(difs .== 1)
            stopmove  = findall(difs .== -1)
        end
        R[k] = p
    end
    C[c] = R
end


# ------------------------------------------------
# PLOT: Corrplot of the variables above
# ------------------------------------------------
# DFr = DFc[rand(1:nrow(DFc),10_000), :]
# @df DFr corrplot([:startWell :stopWell :startWell_tmpl :stopWell_tmpl], 
#     grid = false, method = :pearson, order = :hclust)
