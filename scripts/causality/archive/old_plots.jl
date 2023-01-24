## from plot_alltimes.jl

# ----------------------*
# BIN THE MEANS -- DIFF |
# ----------------------*

B = bin_the_curves(values(Cmd_ca1pfc)...; bins, x_time)
V = collect(skipmissing(values(first(values(Cd_ca1pfc)))))
J = [leaveoneout(V; func=func_bin) for V in collect(values(Cd_ca1pfc))]
J = [vcat(j...) for j in J]
CIs = [(quantile(j[:,col],0.025), quantile(j[:,col], 0.975)) for j in J, col in 1:3]
yerrors = [abs.(c.-b) for (c,b) in zip(CIs,B)]
G = if trans == false
    groupedbar(replace.(collect(keys(Cmd_ca1pfc))," "=>"\n"), B;
               errorbar=yerrors, linewidth=2, label=("early","intermediate","late"),
               legend=:none,
               alpha=0.5, ylabel="Binned 𝔸")
else
    manual_trans(x) = [x[i,j] for j in 1:size(x,2), i in 1:size(x,1)]
    #group = vec(repeat(collect(keys(Cmd_ca1pfc)), outer=(1,3)))
    #nam = repeat(["early","intermediate","late"][Utils.na,:], outer=(4,1))
    #G2 = groupedbar(vec(nam), Matrix(B');
    #           errorbar=[yerrors[i,j] for j in 1:size(yerrors,2), i in 1:size(yerrors,1)], 
    #           linewidth=2, group=vec(group),
    #           alpha=0.5, ylabel="Binned 𝔸")
    #groupedbar(vec(nam), Matrix(B');
    #           errorbar=[yerrors[i,j] for j in 1:size(yerrors,2), i in 1:size(yerrors,1)], 
    #           linewidth=2, group=vec(group), legend=:none,
    #           alpha=0.5, ylabel="Binned 𝔸")
    groupedbar(["early","intermediate","late"], manual_trans(B);
               legend=:none,
               errorbar=yerrors, linewidth=2, 
               grid=false,minorgrid=false,
               label=(string.(collect(keys(Cmd_pfcca1)))),
               alpha=0.5, ylabel="Binned 𝔸")
end
Plot.save("diff - ca1pfc - groupedbar")

B = bin_the_curves(values(Cmd_pfcca1)...; bins, x_time)
V = collect(skipmissing(values(first(values(Cd_pfcca1)))))
J = [leaveoneout(V; func=func_bin) for V in collect(values(Cd_pfcca1))]
J = [vcat(j...) for j in J]
CIs = [(quantile(j[:,col],0.025), quantile(j[:,col], 0.975)) for j in J, col in 1:3]
yerrors = [abs.(c.-b) for (c,b) in zip(CIs,B)]
G = if trans == false
    groupedbar(replace.(collect(keys(Cmd_pfcca1))," "=>"\n"), B;
               legend=:none,
               errorbar=yerrors, linewidth=2, label=("early","intermediate","late"),
               alpha=0.5, ylabel="Binned 𝔸")
else
    manual_trans(x) = [x[i,j] for j in 1:size(x,2), i in 1:size(x,1)]
    #group = vec(manual_trans(repeat(collect(keys(Cmd_pfcca1))[:, Utils.na], outer=(1,3))))
    #nam = repeat(["early","intermediate","late"][Utils.na,:], outer=(4,1))
    #groupedbar(vec(nam), Matrix(B');
    #           errorbar=manual_trans(yerrors), 
    #           #legend=:none,
    #           linewidth=2, group=vec(group),
    #           alpha=0.5, ylabel="Binned 𝔸")
    groupedbar(["early","intermediate","late"], 
               manual_trans(B);
               legend=:none,
               errorbar=yerrors, linewidth=2, 
               label=(string.(collect(keys(Cmd_pfcca1)))),
               alpha=0.5, ylabel="Binned 𝔸")
end
Plot.save("diff - pfcca1 - groupedbar")

