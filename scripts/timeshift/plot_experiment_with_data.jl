


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

Is = @subset(I, :datacut .== datacut, :marginal .== marginal)
Ss = @subset(S, :datacut .== datacut, :marginal .== marginal)
groups = [:datacut, :marginal, :unit, :shift, :area]
G = Table.group.mtg_via_commonmap(groups, Is, Ss)

P_sig, P_non = [], []
Ic = DataFrame()
@showprogress for (ig, sg) in G
    cig = DataFrame(copy(ig))
    # significant
    C = ecdf(sg.value)
    cig[!,:frac] .= mean(cig.value .> sg.value)
    cig[!,:sig]  .= 1 .- cig.frac
    # Plot
    p =plot(C, title="unit=$(Int64(ig.unit[1])) area=$(ig.area[1])",
                        xlabel="bits per spike", ylabel="fraction", ylim=[0,1])
    plot!([ig.value..., ig.value...], [ylims()...])
    if cig.sig[1] < 0.05
        push!(P_sig, p)
    else
        push!(P_non, p)
    end
    # Add to DF
    append!(Ic, cig)
end

using Term
tcolor = hex(colorant"red")

function print_shuf_stats(subset_args...; name="Statistics")
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
end;

bonf_cells_are_sig, cells_are_sig = print_shuf_stats();
print_shuf_stats(:area => x->x.=="CA1"; name="CA1");
print_shuf_stats(:area => x->x.=="PFC"; name="PFC");

print_shuf_stats(:area => x->x.=="CA1", :shift=>x->x.<0; name="CA1 past");
print_shuf_stats(:area => x->x.=="CA1", :shift=>x->x.>0; name="CA1 future");
print_shuf_stats(:area => x->x.=="PFC", :shift=>x->x.>0; name="PFC future");
print_shuf_stats(:area => x->x.=="PFC", :shift=>x->x.<0; name="PFC past");

function plot_info_amount(Is, sig_frac; title="", c=:teal)
    subset_of_interest = subset(bonf_cells_are_sig, :sig => x -> x .== sig_frac)
    insig_cells = Utils.ismember(Is.unit, subset_of_interest.unit)
    #println("Insig cell count  = $(round(sum(insig_cells)/length(unique(Is.shift))))")
    insig_cells = Is[insig_cells, :]
    # Candidate for a Util
    insig_shape = combine(groupby(insig_cells, :shift), 
            :value => mean,
            :value => median,
            :value => (x->nanquantile(x, 0.05) ) => :lower,
            :value => (x->nanquantile(x, 0.95)) => :upper)
    p = @df insig_shape plot(:shift, :value_mean; c,
                              fill=(:lower, :upper), title, label="mean")
    @df insig_shape plot!(:shift, :lower, style=:dash, c=:black, label="lower")
    @df insig_shape plot!(:shift, :upper, style=:dash, c=:black, label="upper")
    p
end

p1 = plot_info_amount(Is, true,  title="Significant cell average")
p2 = plot_info_amount(Is, false, title="Insignificant cell average")
p3 = plot_info_amount(Ss, true,  title="Significant shuffle cell average", c=:red)
p4 = plot_info_amount(Ss, false, title="Insignificant shuffle cell average", c=:red)

P = plot(p1,p2, p3, p4, link=:all)


# -----------
# Corrected
# -----------

function correction(sig_frac) # TODO THINK ABOUT THIS

    sig_of_interest   = subset(bonf_cells_are_sig, :sig => x -> x .== sig_frac)
    insig_of_interest = subset(bonf_cells_are_sig, :sig => x -> x .== sig_frac)

    # TODO THINK ABOUT THIS

    C = Utils.ismember(Is.unit, subset_of_interest.unit)
    C = Is[C, :]

    C = Utils.ismember(Is.unit, subset_of_interest.unit)
    C = Is[C, :]
end
