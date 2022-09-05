
Plot.setfolder("behavior","trajectory-diversity")

s = []
for (g,group) in enumerate(groupby(beh, :stopWell))
    push!(s,scatter(group.x, group.y, label="stopwell=$(group.stopWell[1])", seriescolor=g))
end
plot(s..., markersize=1, markerstrokewidth=0)

sep = []
for group in groupby(beh, :stopWell)
    push!(sep, scatter(group.x, group.y, label="stopwell=$(group.stopWell[1])", group=group.epoch))
end
plot(sep..., markersize=1, markerstrokewidth=0, legend=:none)

sep = []

# Measure angular variability
# ---------------------------
tsk = Load.load_task("RY16",36)
celleps = combine(groupby(spikes, :unit), :epoch => unique, renamecols=false)
Load.filtreg.register(tsk,celleps,on="epoch", transfer=["task"])

# Get occupancy to get angular var
WIDTHS = OrderedDict(
  "x"=>2.5f0, "y"=>2.5f0, 
  "headdir"=>Float32(2pi/80),
  "stopWell"=>1
)

# Measure angular variability
# ---------------------------
widths = OrderedDict(zip(keys(WIDTHS), values(WIDTHS).*1))
props     = ["x","y","headdir","stopWell"]
radiusinc = [0.1f0, 0.1f0, 0.01f0, 0f0]
thresh = 1f0
maxrad = nothing
@time G = Field.adaptive.get_grid(beh, props; widths, thresh, maxrad, radiusinc);
O = @time Field.adaptive.get_occupancy(beh, G);

C = Float32.(copy(O.count))
view(C, repeat(all(C .== 0, dims=4), 1,1,1,7)) .= NaN32

pathvar= []
for c in eachslice(C,dims=4)
    push!(pathvar,nanvar(c,dims=3))
end
for i in 1:length(pathvar)
    pathvar[i] = Utils.squeeze(pathvar[i])
end
H = []
for (stopWell,i) in zip(-1:5, 1:length(pathvar))
    push!(H,heatmap(log10.(pathvar[i]), title="$stopWell"))
end
plot(H...)
Plot.save("Path variance per xy (var(headdir counts)) per stopWell")


bar(-1:5, [nanmedian(p) for p in pathvar], xlabel="StopWell", ylabel="Head direction variance during runs")
Plot.save("Median path variance per xy (var(headdir counts)) per stopWell")

[nanextrema(p) for p in pathvar]



# VARIABILITY IN STATE TRANSITION
# -------------------------------

