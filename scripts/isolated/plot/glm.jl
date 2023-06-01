include("../imports_isolated.jl")
include("../load_cyclewise_checkpoint.jl")

#    _  _     ____  _       _   _   _             
#  _| || |_  |  _ \| | ___ | |_| |_(_)_ __   __ _ 
# |_  ..  _| | |_) | |/ _ \| __| __| | '_ \ / _` |
# |_      _| |  __/| | (_) | |_| |_| | | | | (_| |
#   |_||_|   |_|   |_|\___/ \__|\__|_|_| |_|\__, |
#                                           |___/ 

Plot.setfolder("isolated", "glm")
Plot.setappend((;animal=opt["animal"], day=opt["day"], tet=opt["tet"]))
Plot.printstate()

# modelsort(x) = sort(OrderedDict(x), by= q->[q.unit, q.relcyc])
# It = Iterators.product(zip(modelsort.((model_cellhasiso,)), 
#             ("cell has iso", )), ["adjr2", "mae"])
#
# P = Dict()
# for ((model, name), measure) in It
#     D = to_dataframe(model, grabfield(measure))
#     Dsum = sort(combine(groupby(D, [:relcyc, :indep]), 
#         :value=> (x->nanmean(collect(DIutils.skipnan(x)))) =>:value),
#             [:indep,:relcyc])
#     title = "$name $measure"
#     p = begin
#         @df @subset(Dsum,:indep .== "pfcca1") begin
#             plot(:relcyc, :value, alpha=0.5, label="pfc->ca1")
#         end
#         @df @subset(Dsum,:indep .== "ca1pfc") begin
#             plot!(:relcyc, :value, alpha=0.5, label="ca1->pfc",
#                 ylabel=measure)
#         end
#     end
#     push!(P, title=>p)
# end
# plot(values(P)...)
#
# # View predictions
#
#
# using Blink, Interact, Interpolations
# function interp(D)
#     D = copy(D)
#     x,y,z = 1:size(D,1), D.value, D.relcyc
#     notmissing = (!).(isnan.(y))
#     xn, yn, zn = x[notmissing], y[notmissing], z[notmissing]
#     itp = linear_interpolation(xn, yn)
#     itpr = linear_interpolation(xn, zn)
#     x = unique(clamp.(x, xn[1], xn[end]))
#     D[x,:value] = (itp[x,:])
#     D[x,:relcyc] = (itpr[x,:])
#     D = D[(!).(isnan.(D.value)), :]
#     D
# end
#
# using Blink, Intereact
# for ((model, nm), measure) in It
#
#     tle="$nm $measure"
#     D = sort(to_dataframe(model, 
#                  grabfield(measure)), [:relcyc])
#     P=[(#d=interp(d);
#         @df d (
#         scatter(:relcyc, :value; alpha=0.5, ylabel=measure, 
#         title=tle, linestyle=:solid, tickfontsize=6, label=replace(:indep[1],
#             "pfcca1"=>"pfc ⟶ ca1",
#             "ca1pfc"=>"ca1 ⟶ pfc"),
#     c=:indep[1] .== "pfcca1" ? :red : :blue)); 
#     #@df d plot!(:relcyc, :value)
#     )
#     for d in groupby(D, :unit)]
#     ui = @manipulate for i in eachindex(P)
#         P[i]
#     end
#     w=Window(); body!(w,ui)
#
# end
#
# plot(
#     plot(
#         glmplot(model_isocount, "pfcca1", grabfield, label="iso count"),
#         glmplot(model_hasiso,   "pfcca1", grabfield, label="iso occurance");
#         link=:y,
#         layout=grid(2,1),
#         title= "pfc -> ca1"
#     ),
#     plot(
#         glmplot(model_isocount, "ca1pfc", grabfield, label="iso count"),
#         glmplot(model_hasiso,   "ca1pfc", grabfield, label="iso occurance");
#         link=:y,
#         layout=grid(2,1),
#         title="ca1 -> pfc"
#     ),
#     ylims=(0,0.4),
#     size=(600, 900)
# )
# Plot.save("iso, count and occur, date=$(Dates.now())")
#
#
# # # # # # # ### #### Shuffles #  # #################
# model, shuf = model_cellhasiso, shuffle_cellhasiso
# shuf = Dict((;hash=k)=>v for (k,v) in shuffle_cellhasiso)
# measure = "mae"
# D  = Table.to_dataframe(model,   grabfield(measure))
# Ds = Dict(k => ( Table.to_dataframe(v, grabfield(measure)))
#         for (k,v) in shuf)
# Ds = Table.to_dataframe(Ds)
#
# import Gadfly
# Gadfly.plot(combine(groupby(D, :relcyc), :value=>nanmedian, renamecols=false),  
#     x=:relcyc, y=:value, Gadfly.Geom.point)
# Gadfly.plot(Ds, x=:relcyc, y=:value, Gadfly.Theme(alphas=[0.01]), 
#             Gadfly.Geom.point)
#
#
# histogram(D.value, alpha=0.5, normalize=:pdf)
# histogram!(Ds.value, alpha=0.5, normalize=:pdf)
# vline!([nanmedian(D.value)], c=:blue)
# vline!([nanmedian(Ds.value)], c=:red, linestyle=:dash)
# # xlims!(0.005,0.01)
#
# # using RCall
# # @rlibrary ggplot2
# # ggplot(D, aes(x="relcyc", y="value")) + 
# #     stat_smooth(color="blue")
#
# function show_spiking_std(df::DataFrame; isolated=false)
#     cols = string.(filter( x-> x!==nothing, tryparse.(Int,names(df))))
#     if isolated
#         cols = cols .* "_i"
#         cols = intersect(cols, names(df))
#     end
#     n= df[!,cols] 
#     N = Matrix(n)
#     heatmap(N./std(N,dims=1), clims=(0,20), c=:grays)
# end
# show_spiking_std(df)
# show_spiking_std(df;isolated=true)
#
#
# function describe_cycstats(ca1cycstat, pfccycstat)
#     kws=(;label="")
#     p1=plot(
#      (@df ca1cycstat histogram(:isolated_sum;xlabel="iso spikes emitted",kws...)),
#      (@df ca1cycstat histogram(:isodiv;xlabel="# of uniq iso cells",kws...)),
#      (@df ca1cycstat histogram(:adjdiv;ylabel="# of uniq adj cells",kws...)),
#      Plot.blank((plot();Plots.annotate!(0.5,0.5,text("CA1 cycle\nstatistics",14))),
#             visible=false, size=(100,50)),
#      layout=grid(2,2));
#     p2=plot(
#      (@df pfccycstat histogram(:pfcisosum;xlabel="iso spikes emitted",kws...)),
#      (@df pfccycstat histogram(:pfcisodiv;xlabel="# of uniq iso cells",kws...)),
#      (@df pfccycstat histogram(:pfcadjdiv;ylabel="# of uniq adj cells",kws...)),
#      Plot.blank((plot();Plots.annotate!(0.5,0.5,text("PFC cycle\nstatistics",14))),
#                 visible=false, size=(100,50)),
#         layout=grid(2,2)
#     );
#     plot(p1,p2, size=(1000,500))
# end
# describe_cycstats(ca1cycstat, pfccycstat)

function remove_y_ypred_yr_ypredr_keys(model::Dict)
    @showprogress "removing y keys" for (k,v) in model
        if haskey(v, :y)
            delete!(v, :y)
        end
        if haskey(v, :ypred)
            delete!(v, :ypred)
        end
        if haskey(v, :yr)
            delete!(v, :yr)
        end
        if haskey(v, :ypredr)
            delete!(v, :ypredr)
        end
    end
    return model
end

# ----------------------------
# Ready out data frames
# ----------------------------
    # model   = initorget("model_cellhasiso"; obj=Dict())
    model   = remove_y_ypred_yr_ypredr_keys(model_spikecount)
    shuffle = remove_y_ypred_yr_ypredr_keys(shuffle_spikecount)
    shuffle   = initorget("shuffle_spikecount"; obj=Dict())

    model   = remove_y_ypred_yr_ypredr_keys(model_cellhasiso)
    shuffle = remove_y_ypred_yr_ypredr_keys(shuffle_cellhasiso)
    # shuffle   = initorget("shuffle_spikecount"; obj=Dict())

df_iso  = @time Table.to_dataframe(model, name="value", key_name="prop")
dfs_iso = @time Table.to_dataframe(shuffle, name="value", key_name="prop")

df_iso_val  = @time @subset(df_iso, :prop .== "pseudo_R2")
dfs_iso_val = @time @subset(dfs_iso, :prop .== "pseudo_R2")

commit_cycwise_vars(["df_iso", "df_iso_val", "dfs_iso", "dfs_iso_val"])
# ----------------------------

# ----------------------------
# General distribution of values
# ----------------------------
    # RAW VALUES plot distribution of values
    @df df_iso_val begin
        histogram(:value, bins=100, alpha=0.5, label="pseudo_R2",
            title="positive fraction = $(mean(:value .> 0))")
    end
    @df df_iso_val begin
    ecdfplot!(collect(DIutils.skipnan(:value)), bins=100, alpha=0.5,
        label="pseudo_R2",
        normalize=:pdf, title="positive fraction = $(mean(:value .> 0))")
    end
    vline!([0], c=:red)
# ----------------------------

# ----------------------------
# Relative cycles
# ----------------------------
    # Rel1ative cycle plot distribution of values
    @df df_iso_val begin
        scatter(:relcyc, :value, alpha=0.1, label="pseudo_R2",
            title="positive fraction = $(mean(:value .> 0))")
    end
    @df combine(groupby(df_iso_val, :relcyc), :value=>nanmean, renamecols=false) begin
        scatter!(:relcyc, :value, alpha=0.1, label="pseudo_R2")
    end
    vline!([0.0], c=:red)
    hline!([0.0], c=:red)
# ----------------------------

# ----------------------------
# 10 MINUTE BINS plot distribution of values
# ----------------------------
    # Bin plot distribution of values
    @df @subset(df_iso, :prop .== "pseudo_R2") begin
        scatter(:bin, :value, alpha=0.1, label="pseudo_R2",
            title="positive fraction = $(mean(:value .> 0))")
    end

    # BINS WITH POSITIVE PSEUDO R2
    d = @subset(df_iso, :prop .== "pseudo_R2", :value .> 0.0)
    scatter(d.bin, d.value, group=d.unit, alpha=0.5, label="",
        title="positive fraction = $(mean(:value .> 0))",
        ylims=(-0.001, 0.1))
    t = Plots.text.(d.dir, [8], ["black"])
    p=annotate!.(d.bin, d.value, t)
# ----------------------------


# -----------------------------------------------------
# 10 MINUTE BINS plot distribution of values via ggplot
# -----------------------------------------------------
    using RCall
    @rlibrary ggplot2
    @rimport ggplot2
    @rlibrary ragg

    dfi = @subset(df_iso, :dist .== "binomial")
    begin
        # RELCYC VERSUS PSEUDO R2
        ggplot(dfi, aes(x=:relcyc,y=:value)) + 
               stat_bin(bins=100, geom="line", alpha=0.5) +
               ggplot2.theme(family="serif", text=element_text(size=20))

        # BINS VERSUS PSUEDO R2
        ggplot(dfi, aes(x=:bin,y=:value)) + 
            stat_bin_summary(method = "loess", binwidth = 0.2) +
               # stat_bin(bins=100, geom="line", alpha=0.5) +
               ggplot2.theme(family="Courier", text=element_text(size=8))
    end
# -----------------------------------------------------


# ----------------------------
## DISTRIBUTION BASED
# ----------------------------
P=[]
for dist in unique(df_iso_val.dist)
    relcycs = sort(unique(df_iso_val.relcyc))
push!(P,
        scatter(relcycs, sort(combine(groupby(@subset(df_iso_val,:dist .== dist),
[:relcyc,:dist]), :value=>x->nanmedian(filter(x->x!=-Inf && x > 0, x))),
    :relcyc).value_function, 
-            title=dist, xlabel="relcyc", ylabel="median pseudo_R2")
)
vline!([0.0], c=:red, linestyle=:dash)
end
plot(P..., layout=(length(P), 1))
# ----------------------------


# ----------------------------
# Create shuffle meaned r2 data
# ----------------------------

df_iso_val = @time @subset(df_iso, :prop .== "mae")
dfs_iso_val = @time @subset(dfs_iso, :prop .== "mae")
# Remove nans
df_iso_val = @subset(df_iso_val, :value .!= -Inf .&& :value .!== NaN)
dfs_iso_val = @subset(dfs_iso_val, :value .!= -Inf .&& :value .!== NaN)
mu_dfs_iso_val = combine(groupby(dfs_iso_val, [:unit, :dir, :dist, :relcyc, :bin]), 
    :value=>nanmean, renamecols=false)
H = []
for dir in unique(df_iso_val.dir)
    m, d = @subset(mu_dfs_iso_val, :dir .== dir), 
        @subset(df_iso_val, :dir .== dir)
    push!(H,histogram(m.value, bins=100, alpha=0.5, label="mae shuffle",
        normalize=:probability))
    histogram!(d.value, bins=100, alpha=0.5, label="mae pred",
                normalize=:probability)
end
plot(H..., layout=(length(H), 1))
print(HypothesisTests.KruskalWallisTest(mu_dfs_iso_val.value, df_iso_val.value))

# mean(dfs_iso_val.value), mean(df_iso_val.value)
pred_error = copy(df_iso_val)
G_df  = groupby(pred_error, [:unit, :dir, :dist, :relcyc, :bin])
G_dfs = groupby(mu_dfs_iso_val, [:unit, :dir, :dist, :relcyc, :bin])
keyz = intersect(NamedTuple.(keys(G_df)), NamedTuple.(keys(G_dfs)))
for k in keyz
    @infiltrate
    G_df[k].value = G_dfs[k].value ./ G_df[k].value
end
@subset!(pred_error, :value .!= -Inf, :value .!== NaN, :value .!= Inf)
histogram(pred_error.value, bins=100, alpha=0.5, label="mae prediction gain")

scatter(pred_error.relcyc, pred_error.value, alpha=0.1, 
        label="mae prediction gain",
        title="positive fraction = $(mean(pred_error.value .> 0))")

using NaNStatistics
summ = combine(groupby(pred_error, :relcyc), :value=>
    x->mean(collect(DIutils.skipnan(x))), renamecols=false)
@df summ begin
    plot(:relcyc, :value, alpha=1, label="mae prediction gain",c=:red,
        ylims=(0, 3.0), title="model positive fraction = $(round(mean(pred_error.value .> 0), digits=3))",
    fillrange=0.0, fillalpha=0.1, fillcolor=:red, 
    xlabel="relcyc", ylabel="mean\nmean prediction gain"
    )
end
vline!([0.0], c=:red, linestyle=:dash,label="")
