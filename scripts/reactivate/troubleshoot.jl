include("imports.jl")
include("prep2plot.jl")
using DataFrames

# This analysis indicates that the :component column is dropping templates
H2 = []
# global DFc = subset(DF,  subs[1], view=true)
@time tmpc=combine(groupby(DFc, [:areas, :traj]),
[:startWell_tmpl, :stopWell_tmpl] => 
    ((x,y) -> length(unique(eachrow([x y])))) => :tmpl_combos)
h2=histogram(tmpc.tmpl_combos, group=tmpc.areas, bins=1:1:20, normed=true, 
    alpha=0.2, bar_position=:stack,
    title="histogram of template combos per trajectory", 
    xlabel="number of template combos", ylabel="frequency")
push!(H2, h2)
@showprogress for s in 2:length(subs)
    global DFc = subset(DFc, subs[s], view=true)
    @time tmpc=combine(groupby(DFc, [:areas, :traj]),
        [:startWell_tmpl, :stopWell_tmpl] => 
            ((x,y) -> length(unique(eachrow([x y])))) => :tmpl_combos)
    h2=histogram(tmpc.tmpl_combos, group=tmpc.areas, bins=1:1:20, normed=true, 
        alpha=0.2, bar_position=:stack,
        title="histogram of template combos per trajectory", 
                xlabel="number of template combos", ylabel="frequency")
    push!(H2, h2)
end


# This analysis indicates IF the :component column is equally numbered across
# the templates
combine(groupby(DF, [:startWell_tmpl, :stopWell_tmpl]),
    :component => (x->length(unique(x))) => :count,
    :component => (x->maximum(unique(x))) => :max,
    :component => (x->minimum(unique(x))) => :min
)
tmp=ans tmp.startstop = string.(tmp.startWell_tmpl, "-", tmp.stopWell_tmpl)

b=@df tmp bar(:startstop,   :count, title="Number of components per template")

s=@df tmp scatter(:startstop, :max, title="Number of components per template")
@df tmp scatter!( :startstop,  :min, title="Number of components per template")
# Draw a line between the points
p = plot!()
plot!([p.series_list[1][:x] p.series_list[1][:x]]', [tmp.min tmp.max]', 
    title="Number of components per template", label="", color=:black)



