using DrWatson
using Revise
@time include(scriptsdir("timeshift", "Initialize_plot_experimentation.jl"))
using StatsBase
using Term

using REPL.TerminalMenus
function varrequest(var::Symbol)
    msg=  "Pick a $(string(var))"
    var = unique(S[!,var])
    var[request(msg, RadioMenu(string.(var)))]
end

# Datacut
datacut  = varrequest(:datacut)

# Marginal
marginal = varrequest(:marginal)

# Grab the proper subsets
Is = @subset(I, :datacut .== datacut, :marginal .== marginal)
Ss = @subset(S, :datacut .== datacut, :marginal .== marginal)
groups = [:datacut, :marginal, :unit, :shift, :area]

# Score the dataframe with significance
Is = Timeshift.operation.score_significance(Is, Ss)

function print_shuf_stats(subset_args...; name="Statistics")
    tcolor = hex(colorant"red")
    if !(isempty(subset_args))
        Ics = subset(Ic, subset_args...)
    else
        Ics = Ic
    end
    println(hLine(name))
    cells_are_sig = combine(groupby(Ics, :unit), :sig => (x-> any(x.<0.05)) => :sig )
    frac = mean(cells_are_sig.sig);
    tprintln("Fraction of [#$tcolor]significant[/#$tcolor] cells\n=> [#$tcolor]$(frac)[/#$tcolor]")
    println()
    bonf_cells_are_sig = combine(groupby(Ics, :unit), :sig => (x-> any(x.<0.05/size(Ics,1))) => :sig);
    frac = mean(cells_are_sig.sig);
    tprintln("Bonf. fraction of [#$tcolor]significant[/#$tcolor] cells\n=> [#$tcolor]$(frac)[/#$tcolor]")
    print(hLine(""))
    return bonf_cells_are_sig, cells_are_sig
end

bonf_cells_are_sig, cells_are_sig = print_shuf_stats();
print_shuf_stats(:area => x->x.=="CA1"; name="CA1");
print_shuf_stats(:area => x->x.=="PFC"; name="PFC");

print_shuf_stats(:area => x->x.=="CA1", :shift=>x->x.<0; name="CA1 past");
print_shuf_stats(:area => x->x.=="CA1", :shift=>x->x.>0; name="CA1 future");
print_shuf_stats(:area => x->x.=="PFC", :shift=>x->x.>0; name="PFC future");
print_shuf_stats(:area => x->x.=="PFC", :shift=>x->x.<0; name="PFC past");


p1 = Plot.timeshift.plot_shift_versus_info(Is, true, bonf_cells_are_sig,  title="Significant cell average")
p2 = Plot.timeshift.plot_shift_versus_info(Is, false,bonf_cells_are_sig, title="Insignificant cell average")
p3 = Plot.timeshift.plot_shift_versus_info(Ss, true, bonf_cells_are_sig, title="Significant shuffle\ncell average",   c=:red)
p4 = Plot.timeshift.plot_shift_versus_info(Ss, false,bonf_cells_are_sig, title="Insignificant shuffle\ncell average", c=:red)

P = plot(p1, p2, p3, p4, link=:all)

# -----------
# Fields that are sig
# -----------
# A. Most, zero, and least
subset_of_interest = subset(bonf_cells_are_sig, :sig => x -> x .== true)
P_lmz = Array{Plots.Plot, 2}(undef, size(subset_of_interest,1), 3)
for (u, unit) ∈ enumerate(subset_of_interest.unit)

    ic = @subset(Ic, :unit .== unit)
    ic = sort(ic, :frac)
    fr = round(@subset(cells, :unit .==unit).meanrate[1], digits=1)

    # least sig shift
    #
    # G 
    sₗ = ic[argmin(ic.frac), :shift]
    # most sig shift
    sₘ = ic[argmax(ic.frac), :shift]


p1 = Plot.timeshift.plot_shift_versus_info(Is, true, bonf_cells_are_sig,  title="Significant cell average")
p2 = Plot.timeshift.plot_shift_versus_info(Is, false, title="Insignificant cell average")
p3 = Plot.timeshift.plot_shift_versus_info(Ss, true,  title="Significant shuffle\ncell average",   c=:red)
p4 = Plot.timeshift.plot_shift_versus_info(Ss, false, title="Insignificant shuffle\ncell average", c=:red)

groupby(Is, :unit)
P = plot(p1,p2, p3, p4, link=:all)

# -----------
# Fields that are sig
# -----------
# A. Most, zero, and least
subset_of_interest = subset(bonf_cells_are_sig, :sig => x -> x .== true)
P_lmz = Array{Plots.Plot, 2}(undef, size(subset_of_interest,1), 3)
for (u, unit) ∈ enumerate(subset_of_interest.unit)

    ic = @subset(Ic, :unit .== unit)
    least = @subset(F, :unit .== unit, :shift .== sₗ, :datacut .== datacut, :marginal .== marginal)
    leastᵥ = least.value[1]
    leastₛ = round(least.shift[1], sigdigits=2)

    most = @subset(F, :unit .== unit, :shift .== sₘ, :datacut .== datacut, :marginal .== marginal)

    mostᵥ = most.value[1]
    mostₛ = round(most.shift[1], sigdigits=2)

    zilch = @subset(F, :unit .== unit, :shift .== 0, :datacut .== datacut, :marginal .== marginal)
    zilchᵥ = zilch.value[1]

    clim = cat(leastᵥ, mostᵥ, zilchᵥ, dims=3)
    clim[clim.==0] .= NaN
    clim = ([nanquantile(vec(clim), 0.05), nanquantile(vec(clim), 0.9)]...,)
    #clim=(-Inf,Inf)

    P_lmz[u, 1] = heatmap(leastᵥ; clim, title="area=$(ic.area[1]),\nunit=$unit, fr=$fr, \nleast=$leastₛ")
    P_lmz[u, 2] = heatmap(mostᵥ;  clim, title="area=$(ic.area[1]),\nunit=$unit, fr=$fr, \nmost=$mostₛ")
    P_lmz[u, 3] = heatmap(zilchᵥ; clim, title="area=$(ic.area[1]),\nunit=$unit, fr=$fr, \nzero")

end
w= Window()
ui = @manipulate for u in 1:size(P_lmz, 1)
    plot(P_lmz[u,:]...; layout=grid(1,3), aspect_ratio=1, size=(1000,1000))
end
body!(w, ui)


# B. Most, zero, and least
# All, with most zero and least marked
subset_of_interest = subset(bonf_cells_are_sig, :sig => x -> x .== true)
for cell ∈ subset_of_interest.unit
    # Shifts versus 
    ic = @subset(Ic, :unit .== unit)
    ic = sort(ic, :frac)
    fr = round(@subset(cells, :unit .==unit).meanrate[1], digits=1)
    sort!(ic, :shift)

    bps = @df plot(:shift, :value, ylabel="Bits per spike")
    f   = sort(@subset(F, :datacut .== datacut, 
                       :marginal .== marginal, 
                       :unit .== unit), :shift)

    layout = grid(1, size(f,1))

end

w=Window()
clim = cat([r.value for r in eachrow(f)]..., dims=3)
clim = ([nanquantile(vec(clim), 0.05), nanquantile(vec(clim), 0.9)]...,)
# least sig shift
sₗ = ic[argmin(ic.frac), :]
# most sig shift
sₘ = ic[argmax(ic.frac), :]
ui = @manipulate for r in 1:size(f,1)
    sᵣ = @subset(ic, :shift .== f.shift[r])
    bgcolor, addition = nothing, nothing
    if f[r,:shift] == sₗ.shift
        bgcolor = :red
        i = round(sᵣ.value[1],     sigdigits=2)
        δ = round(i - sₘ.value[1], sigdigits=2)
        addition = "i=$i, δ=$δ"
    elseif f[r,:shift] == sₘ.shift
        i = round(sᵣ.value[1], sigdigits=2)
        δ = round(i - sₗ.value[1], sigdigits=2)
        addition = "i=$i, δ=$δ"
        bgcolor = :green
    else
        bgcolor = :white
        i = round(sᵣ.value[1], sigdigits=2)
        addition = "i=$i"
    end
    bps = @df plot(:shift, :value, ylabel="Bits per spike")
    hm = heatmap(f[r,:].value; clim, title="shift=$(round(f[r,:shift]*60,sigdigits=2))\n$addition", bgcolor)
end
body!(w,ui)


# -----------
# Corrected
# -----------
# very suss effect on the cells. should this correction be cell wise or a
# blanket average corretion? it seems to take over the effect in most cases.
Isc = Timeshift.operation.correct_signalwithshuf(Is, Ss, bonf_cells_are_sig).Icorr
Plot.timeshift.unitshift_heatmap(Isc)

Isc = Timeshift.norm_information(Isc)
Plot.timeshift.unitshift_heatmap(Isc)


