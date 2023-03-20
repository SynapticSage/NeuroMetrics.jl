"""
umap_project.jl

Project the data onto the behavior variables from UMAP -> PCA -> Linear project
onto behavior
"""

Plot.setparentfolder("manifold", "umap_project")
Plot.setappend((;animal, day, filt, s)
Plot.deleteplotfiles()
Plot.printstate()

# beh_vars = [:x, :y, :cuemem, :stopWell, :startWell, :headdir]
beh_vars = [:x, :y, :cuemem]
B = Matrix(beh[:, beh_vars])

# Create a dummy code for behavior stopWells
DS = DIutils.statistic.dummycode(beh[:,:stopWell])
stopWellVars = [Symbol("stopWell_$i") for i in 1:size(DS,2)]
B = hcat(B, DS)
# Create a dummy code for behavior startWells
DS = DIutils.statistic.dummycode(beh[:,:startWell])
startWellVars = [Symbol("startWell_$i") for i in 1:size(DS,2)]
B = hcat(B, DS)

# Subset for this data fraction
B = B[inds_of_t[s], :]

C(x;subset=Colon(),kws...)=DIutils.plotutils.getplotcolor(beh[subset,x], :vik; kws...)
c=DIutils.plotutils.getplotcolor(beh[:,:cuemem], :vik)

# """
#         project_onto_behavior(R::Matrix, beh::DataFrame, beh_vars::Vector)
# Project the data `R` onto the behavior variables `beh_vars` in `beh`.
# # Arguments
# - `R::Matrix`: The data to project, time x neurons
# - `beh::DataFrame`: The behavior data, time x properties
# # Returns
# - `Rproj::Matrix`: The projected behavior data, time x properties
# """
# function project_onto_behavior(r::AbstractMatrix, beh::DataFrame, times, beh_vars::Vector)
#     B = Matrix(beh[times, beh_vars])
#     # Matrix(beh[times, beh_vars]) * r' * inv(r * r') * r
#     r * B' * pinv(B * B') * B
# end

# PCA both into same dims
pca = fit(PCA, Float64.(em), maxoutdim=3, pratio=1.0);
pca2 = fit(PCA, Float64.(Matrix(B)), maxoutdim=3, pratio=1.0);

# Normalize columsn to 0 to 1
# Center
function SC(B)
    B = (B .- minimum(B, dims=1)) ./ (maximum(B, dims=1) .- minimum(B, dims=1))
    # Center 
    B = B .- mean(B, dims=1)
end
pca2all = fit(PCA, SC(B), maxoutdim=3, pratio=1.0);

plot(eachcol(pca2all.proj)...;c, alpha=0.25, legend=false, 
    xlabel="PC1", ylabel="PC2", zlabel="PC3", title="Behavior")
Plot.save("behavior_pca_colorby_cuemem")
plot(eachcol(pca.proj)...; c,
    xlabel="PC1", ylabel="PC2", zlabel="PC3", title="EM")
Plot.save("em_pca_colorby_cuemem")


using MultivariateStats
lsq = llsq(Float64.(em), B[inds_of_t[s], :])
emp = (em * lsq[1:end-1,:]) .+ DIutils.arr.atleast2d(lsq[end, :])'
lsq_coef = lsq[1:end-1,:]
lsq_bias = lsq[end, :]
pca_emp = fit(PCA, emp, maxoutdim=3, pratio=1.0);
plot(
plot(eachcol(pca_emp.proj)...; c),
scatter(beh[:,:cuemem]; c, alpha=0.1)
)


mani_kws = (; markersize=0.5, markerboderstyle=:dot, markerstrokealpha=0.01,
    markerstrokewidth=0.1,)
scatter_kws = (; markersize=3, markerboderstyle=:dot, markerstrokealpha=0.01,
    markerstrokewidth=0.1,)
beh[:, :speed_clip] = clamp.(beh[:, :speed], 0, 5)
B = beh[inds_of_t[s], :]

x = :speed_clip
plot(
    plot(eachcol(pca_emp.proj)...; c=C(x;alpha=:lowtohigh)),
    scatter(beh[:,x]; c=C(x), alpha=0.1),
    title=string(x)
)
Plot.save("Behaviorally relavant, color by $(x)")

x = :blocktraj
plot(
    plot(eachcol(pca_emp.proj)...; c=C(x;alpha=:none)),
    scatter(beh[:,x]; c=C(x), alpha=0.1),
    title=string(x)
)
Plot.save("Behaviorally relavant, color by $(x)")

x = :blocktraj
beh_sub = findall(B[!, :cuemem] .!= -1)
plot(
    plot(eachcol(pca_emp.proj[beh_sub,:])...; c=C(x;subset=beh_sub,alpha=:none)),
    scatter(B[beh_sub,x]; c=C(x,subset=beh_sub,), alpha=0.1),
    title=string(x)
)
Plot.save("Behaviorally relavant, color by $(x)")

# Get dict translating haatrajnum to hatraj
haatrajnum_to_hatraj = Dict(
    k=>v for (k,v) in zip(
        unique(beh[:, :hatrajnum]),
        unique(beh[:, :hatraj])
    )
    if !ismissing(k) && !ismissing(v)
)

x = :hatrajnum
B = beh[inds_of_t[s], :]
beh_sub = findall((!).(ismissing.(B[!, :hatrajnum])))
plot(
    plot(eachcol(pca_emp.proj[beh_sub,:])...; mani_kws...,
        c=C(x;subset=beh_sub,alpha=:none)),
    # Set xticks to hatraj string labels
    scatter(B[beh_sub,x]; c=C(x,subset=beh_sub,), alpha=0.1,
            scatter_kws...,
            yticks=(
                collect(keys(haatrajnum_to_hatraj)),
        collect(values(haatrajnum_to_hatraj))),
    ),
    title=string(x)
)
# Plot.save("Behaviorally relavant, color by $(x)")

# Color by stopWell
x = :stopWell
B = beh[inds_of_t[s], :]
beh_sub = findall((!).(
    ismissing.(B[!, :stopWell]) .&&
    B.stopWell .!= -1
))
plot(
    scatter(eachcol(pca_emp.proj[beh_sub,:])...; 
        c=C(x;subset=beh_sub,alpha=:none),
        mani_kws...
    ),
    # Set xticks to hatraj string labels
    scatter(B[beh_sub,x]; c=C(x,subset=beh_sub,), alpha=0.1,
        scatter_kws...
    ),
    title=string(x)
)

# Color by startWell
x = :startWell
B = beh[inds_of_t[s], :]
beh_sub = findall((!).(
    ismissing.(B[!, :startWell]) .&&
    B.startWell .!= -1
))
plot(
    scatter(eachcol(pca_emp.proj[beh_sub,:])...; 
        c=C(x;subset=beh_sub,alpha=:none),
        mani_kws...
    ),
    # Set xticks to hatraj string labels
    scatter(B[beh_sub,x]; c=C(x,subset=beh_sub,), alpha=0.1,
        scatter_kws...
    ),
    title=string(x),
    background_color=:grey90
)

# Create a new field enumerating (startWell, stopWell) combinations
# with integers
S = collect(eachrow(Matrix(beh[:, [:startWell, :stopWell]])))
uS = unique(S, dims=1)
sf = [findfirst(isequal(s), uS) for s in S]
beh[:, :startstopWell] = sf

B = beh[inds_of_t[s], :]
x = :startstopWell
beh_sub = findall((!).(
    ismissing.(B[!, x]) .&&
    B[!, x] .!= -1
))
plot(
    scatter(eachcol(pca_emp.proj[beh_sub,:])...; 
        c=C(x;subset=beh_sub,alpha=:none),
        mani_kws...
    ),
    # Set xticks to hatraj string labels
    scatter(B[beh_sub,x]; c=C(x,subset=beh_sub,), alpha=0.1,
        scatter_kws...
    ),
    title=string(x),
    background_color=:grey90
)

# julia> cor([emp Matrix(beh[inds_of_t[s], [:x,:y,:cuemem]]) ])
# 6×6 Matrix{Float64}:
#   1.0         -0.0399363   0.966072     0.16858     -0.00230479   0.257325
#  -0.0399363    1.0         0.146473    -0.00673286   0.0577124    0.0390151
#   0.966072     0.146473    1.0          0.16286      0.00845338   0.266362
#   0.16858     -0.00673286  0.16286      1.0         -0.333641     0.147966
#  -0.00230479   0.0577124   0.00845338  -0.333641     1.0         -0.161629
#   0.257325     0.0390151   0.266362     0.147966    -0.161629     1.0
# Correlation with behavior variables not amazing

# Get color of cuemem
# Just looking at the correlation between the behavior variables 
for beh_var in beh_vars
    println("beh_var=$beh_var")
    println(begin
        combine(groupby(
            DataFrame([beh[inds_of_t[s], beh_var] emp[:,:]], :auto), :x1),
        [:x2,:x3,:x4] .=> mean)
    end
)
end

[emp B]


# Cannonical corr
out = fit(CCA, Float64.(pca.proj), pca2.proj)

