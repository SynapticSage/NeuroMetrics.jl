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

modelsort(x) = sort(OrderedDict(x), by= q->[q.unit, q.relcyc])
It = Iterators.product(zip(modelsort.((model_cellhasiso,)), 
            ("cell has iso", )), ["adjr2", "mae"])

P = Dict()
for ((model, name), measure) in It
    D = to_dataframe(model, grabfield(measure))
    Dsum = sort(combine(groupby(D, [:relcyc, :indep]), 
        :value=> (x->nanmean(collect(DIutils.skipnan(x)))) =>:value),
            [:indep,:relcyc])
    title = "$name $measure"
    p = begin
        @df @subset(Dsum,:indep .== "pfcca1") begin
            plot(:relcyc, :value, alpha=0.5, label="pfc->ca1")
        end
        @df @subset(Dsum,:indep .== "ca1pfc") begin
            plot!(:relcyc, :value, alpha=0.5, label="ca1->pfc",
                ylabel=measure)
        end
    end
    push!(P, title=>p)
end
plot(values(P)...)

# View predictions


using Blink, Interact, Interpolations
function interp(D)
    D = copy(D)
    x,y,z = 1:size(D,1), D.value, D.relcyc
    notmissing = (!).(isnan.(y))
    xn, yn, zn = x[notmissing], y[notmissing], z[notmissing]
    itp = linear_interpolation(xn, yn)
    itpr = linear_interpolation(xn, zn)
    x = unique(clamp.(x, xn[1], xn[end]))
    D[x,:value] = (itp[x,:])
    D[x,:relcyc] = (itpr[x,:])
    D = D[(!).(isnan.(D.value)), :]
    D
end

for ((model, nm), measure) in It

    tle="$nm $measure"
    D = sort(to_dataframe(model, 
                 grabfield(measure)), [:relcyc])
    P=[(#d=interp(d);
        @df d (
        scatter(:relcyc, :value; alpha=0.5, ylabel=measure, 
        title=tle, linestyle=:solid, tickfontsize=6, label=replace(:indep[1],
            "pfcca1"=>"pfc ⟶ ca1",
            "ca1pfc"=>"ca1 ⟶ pfc"),
    c=:indep[1] .== "pfcca1" ? :red : :blue)); 
    #@df d plot!(:relcyc, :value)
    )
    for d in groupby(D, :unit)]
    ui = @manipulate for i in eachindex(P)
        P[i]
    end
    w=Window(); body!(w,ui)

end

plot(
    plot(
        glmplot(model_isocount, "pfcca1", grabfield, label="iso count"),
        glmplot(model_hasiso,   "pfcca1", grabfield, label="iso occurance");
        link=:y,
        layout=grid(2,1),
        title= "pfc -> ca1"
    ),
    plot(
        glmplot(model_isocount, "ca1pfc", grabfield, label="iso count"),
        glmplot(model_hasiso,   "ca1pfc", grabfield, label="iso occurance");
        link=:y,
        layout=grid(2,1),
        title="ca1 -> pfc"
    ),
    ylims=(0,0.4),
    size=(600, 900)
)
Plot.save("iso, count and occur, date=$(Dates.now())")


# # # # # # ### #### Shuffles #  # #################
model, shuf = model_cellhasiso, shuffle_cellhasiso
shuf = Dict((;hash=k)=>v for (k,v) in shuffle_cellhasiso)
measure = "mae"
D  = Table.to_dataframe(model,   grabfield(measure))
Ds = Dict(k => ( Table.to_dataframe(v, grabfield(measure)))
        for (k,v) in shuf)
Ds = Table.to_dataframe(Ds)

import Gadfly
Gadfly.plot(combine(groupby(D, :relcyc), :value=>nanmedian, renamecols=false),  
    x=:relcyc, y=:value, Gadfly.Geom.point)
Gadfly.plot(Ds, x=:relcyc, y=:value, Gadfly.Theme(alphas=[0.01]), 
            Gadfly.Geom.point)


histogram(D.value, alpha=0.5, normalize=:pdf)
histogram!(Ds.value, alpha=0.5, normalize=:pdf)
vline!([nanmedian(D.value)], c=:blue)
vline!([nanmedian(Ds.value)], c=:red, linestyle=:dash)
# xlims!(0.005,0.01)

# using RCall
# @rlibrary ggplot2
# ggplot(D, aes(x="relcyc", y="value")) + 
#     stat_smooth(color="blue")

function show_spiking_std(df::DataFrame; isolated=false)
    cols = string.(filter( x-> x!==nothing, tryparse.(Int,names(df))))
    if isolated
        cols = cols .* "_i"
        cols = intersect(cols, names(df))
    end
    n= df[!,cols] 
    N = Matrix(n)
    heatmap(N./std(N,dims=1), clims=(0,20), c=:grays)
end
show_spiking_std(df)
show_spiking_std(df;isolated=true)


function describe_cycstats(ca1cycstat, pfccycstat)
    kws=(;label="")
    p1=plot(
     (@df ca1cycstat histogram(:isolated_sum;xlabel="iso spikes emitted",kws...)),
     (@df ca1cycstat histogram(:isodiv;xlabel="# of uniq iso cells",kws...)),
     (@df ca1cycstat histogram(:adjdiv;ylabel="# of uniq adj cells",kws...)),
     Plot.blank((plot();Plots.annotate!(0.5,0.5,text("CA1 cycle\nstatistics",14))),
            visible=false, size=(100,50)),
     layout=grid(2,2));
    p2=plot(
     (@df pfccycstat histogram(:pfcisosum;xlabel="iso spikes emitted",kws...)),
     (@df pfccycstat histogram(:pfcisodiv;xlabel="# of uniq iso cells",kws...)),
     (@df pfccycstat histogram(:pfcadjdiv;ylabel="# of uniq adj cells",kws...)),
     Plot.blank((plot();Plots.annotate!(0.5,0.5,text("PFC cycle\nstatistics",14))),
                visible=false, size=(100,50)),
        layout=grid(2,2)
    );
    plot(p1,p2, size=(1000,500))
end
describe_cycstats(ca1cycstat, pfccycstat)
