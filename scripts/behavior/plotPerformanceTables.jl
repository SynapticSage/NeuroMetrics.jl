using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))

using Measures
import MAT, MATLAB
using DataFrames, CSV
using ProgressMeter
using Plots
import Plots: plot
using Gadfly, Cairo, Fontconfig
includet("/home/ryoung/Code/projects/goal-code/src/utils/SearchSortedNearest.jl/src/SearchSortedNearest.jl")
table_folder = "/home/ryoung/Projects/goal-code/data/exp_raw/performance"
blacklist=Dict(
   "animal"  => ["RY20", "RY22", "RY16", "RY16", "RY16", "RY22"],
     "day"   => [5,      23,     50,     52,     55,     13]
    )
blacklist=DataFrame(blacklist)
mkdir(plotsdir("behavior","implant_performance"))


normmin(time) = time./60 .- minimum(time)/60;

function nanfill(x::AbstractVector, tq; mode="f")
    if all(isnan.(x))
        xq=x
    else
        inds = (!).(isnan.(x));
        x = x[inds]
        t = tq[inds]
        if mode == "f"
            finds = SearchSortedNearest.searchsortedprevious.([t], tq)
        else
            finds = SearchSortedNearest.searchsortednext.([t], tq)
        end
        xq = x[finds]
    end
    return xq
end
function nanfill!(groups::GroupedDataFrame, cols=Vector{String})
    for (g,group) in enumerate(groups)
        for (col,mod) in Iterators.product(cols,("","_l","_h"))
            #println((g,col))
            col = col*mod
            group[!,col] = nanfill(Vector(group[!,col]), Vector(group[!,"time"]),
                                   mode="f");
        end
    end
    return groups
end


function load_tables(table_folder::String; blacklist=nothing)
    csvfiles = readdir(table_folder)
    D = Vector{DataFrame}([])
    @showprogress for file in csvfiles
        d = CSV.read(joinpath(table_folder,file), DataFrame)
        d[!,"id"] = 1:size(d,1)
        if blacklist != nothing &&
            any(all(reshape(Array(d[1,[:animal,:day]]),1,2) 
                    .== Array(blacklist),dims=2))
            print("Skipping")
            continue
        end
        push!(D,d)
    end
    frame = DataFrame()
    for d in D
        if size(d,2)==0
            continue
        end
        frame = vcat(frame, d, cols=:union)
    end
    return frame
end

function remove_interepoch_time(frame)
    groups = groupby(frame, [:animal,:day,:epoch], sort=true)
    df=DataFrame
    for (g, group) in enumerate(groups)
       if g == 1
           continue
       end
       if all(Array(groups[g][1,[:animal,:day]]) != Array(groups[g-1][1,[:animal,:day]]))
           continue
       end
       println("Subtracting")
       println(
               vcat(df(groups[g][1,[:animal,:day,:epoch]]),
                    df(groups[g-1][1,[:animal,:day,:epoch]]))
              ) 
       println("")
       @assert groups[g].epoch[1] > groups[g-1].epoch[1]
       groups[g][!,"time"] .-= minimum(groups[g][!,"time"]) - maximum(groups[g-1][!,"time"])
       groups[g][!,"id"] .-= minimum(groups[g][!,"id"]) - maximum(groups[g-1][!,"id"])
       group = groups[g]
    end
    return combine(groups, identity)
end

function add_statespace(frame)
    function addresult(group, pair; inds=!)
        result, name = pair
        group[!, name]  .= NaN
        group[!, name*"_l"]  .= NaN
        group[!, name*"_h"]  .= NaN
        group[inds, name]  = result[2:end,1]
        group[inds, name*"_l"] = result[2:end,2]
        group[inds, name*"_h"] = result[2:end,3]
    end

    frame = sort(frame, [:animal, :day, :epoch, :id]);

    # Here, we use our matlab statespace function to capture that for each
    # day group per animal
    groups = groupby(frame, [:animal,:day])
    @showprogress for (g, group) = enumerate(groups)
        println("Group = $g")

        perf = copy(group.correct)
        state = mat"getestprobcorrect( $perf, 0.25, 1)"
        addresult(group, state=>"state")

        inds = group.wm.==1
        wmstate = copy(group.correct)[inds]
        mat"""
        try
            $wmstate = getestprobcorrect( $wmstate, 0.25, 1)
        catch
            $wmstate = nan(size($wmstate,1)+1,3);
        end
        """
        addresult(group, wmstate=>"wmstate", inds=inds)

        inds=group.wm.==0
        cuestate = copy(group.correct)[inds]
        mat"""
        try
            $cuestate = getestprobcorrect( $cuestate, 0.25, 1)
        catch
            $cuestate = nan(size($cuestate,1)+1,3);
        end
        """
        addresult(group, cuestate=>"cuestate", inds=inds)
    end
    mat"close all"
    return groups
end

function normalizetimes!(groups)
    for group in groups
        group[!,"time"] = normmin(group.time)
    end
end

function plot_each_statespace(group)

end

function plot_each_statespace_per_day(group)
end

function overlay_choice_of_statespaces()
end

function plot_statespace(frame, name::String; fillna=true, label=nothing,
        title="", inplace=false, time=nothing, kws...)
    frame = copy(frame)
    if time == nothing
        time = :time
    end
    animal, day = frame[1, ["animal", "day"]]
    state = frame[!, [name, name*"_l", name*"_h"]]
    print(state)
    rib = collect(eachcol(state[:,2:3]))
    rib[1] .-= state[:,1]
    rib[1] = abs.(rib[1])
    rib[2] .-= state[:,1]
    rib = Tuple(rib)
    if label == nothing 
        label = "$animal $day $name"
    end
    if title == nothing 
        title = "$animal $day $name"
    end
    if !inplace
        plot(frame[!,time], state[:,1]; ribbon=rib, label=label, fillalpha=0.2, title=title, kws...)
        hline!([0.8], line=(1mm,:dash,:black), label="")
    else
        plot!(frame[!,time], state[:,1]; ribbon=rib, label=label, fillalpha=0.2, title=title, kws...)
    end
end
function plot_statespace!(frame, name::String ;kws...)
    plot_statespace(frame, name::String ;kws..., inplace=true)
end
function plot_statespace(frame, name::Vector{String}; c=nothing, kws...)
    if c == nothing
        p = plot_statespace(frame, name[1]; kws...)
        for i = 2:length(name)
            plot_statespace!(frame, name[i]; kws...)
        end
    else
        p = plot_statespace(frame, name[1]; c=c[1], kws...)
        for i = 2:length(name)
            plot_statespace!(frame, name[i]; c=c[i], kws...)
        end
    end
    return p
end


frames = load_tables(table_folder; blacklist=blacklist)
frame = sort(frame, [:animal, :day, :epoch, :id]);
plot(frames.time, label="original")
frames = remove_interepoch_time(frames)
plot!(frames.time, label="correction")
groups = add_statespace(frames);
normalizetimes!(groups)

tab = combine(groups, identity)
ftab = combine(nanfill!(groupby(copy(tab),[:animal,:day]), ["state","cuestate","wmstate"]),
               identity)
fgroups = groupby(ftab, [:animal, :day])

sps = ["state","cuestate","wmstate"]
#                                    
# ,---.|         |              o|    
# |---'|    ,---.|--- ,---.     .|    
# |    |    |   ||    `---.     ||    
# `    `---'`---'`---'`---'o    |`---'
#                           `---'     
function peranimal(groups, item; kws...)
    p = (plot_statespace(groups[i], item; fillna=true, time=:id, kws...)
               for i in 1:length(groups))
    p = Plots.plot(p..., margins=0cm, framestyle=:grid)
end

theme(:sand)
pal = theme_palette(:auto)

function peranspec(G)
    overall = peranimal(G, "state",    label="", c=pal[1], title=nothing, xtick=[], ytick=[0, 1], titlefont=font(3), ytickfontsize=4)
    cuestate= peranimal(G, "cuestate", label="", c=pal[2], title=nothing, xtick=[], ytick=[0, 1], titlefont=font(3), ytickfontsize=4)
    memstate= peranimal(G, "wmstate",  label="", c=pal[3], title=nothing, xtick=[], ytick=[0, 1], titlefont=font(3), ytickfontsize=4)
    layout = @layout([a;b c])
    F = Plots.plot(overall, cuestate, memstate, layout=layout)
    overlay = peranimal(G, ["cuestate","wmstate"], c=pal[2:3], title=nothing, label="", titlefont=font(3), ytickfontsize=3, xtick=[])
    return F, overlay
end

F, overlay =peranspec(fgroups)
F
overlay
F.attr[:dpi]=250
savefig(F, plotsdir("behavior","implant_performance","fstatespaces.pdf"))
savefig(F, plotsdir("behavior","implant_performance","fstatespaces.svg"))
savefig(F, plotsdir("behavior","implant_performance","fstatespaces.png"))
overlay.attr[:dpi]=250
savefig(overlay, plotsdir("behavior","implant_performance","fstatespaces_cuememOverlay.pdf"))
savefig(overlay, plotsdir("behavior","implant_performance","fstatespaces_cuememOverlay.svg"))
savefig(overlay, plotsdir("behavior","implant_performance","fstatespaces_cuememOverlay.png"))

F, overlay =peranspec(groups)
F
overlay
F.attr[:dpi]=250
savefig(F, plotsdir("behavior","implant_performance","statespaces.pdf"))
savefig(F, plotsdir("behavior","implant_performance","statespaces.svg"))
savefig(F, plotsdir("behavior","implant_performance","statespaces.png"))
overlay.attr[:dpi]=250
savefig(overlay, plotsdir("behavior","implant_performance","statespaces_cuememOverlay.pdf"))
savefig(overlay, plotsdir("behavior","implant_performance","statespaces_cuememOverlay.svg"))
savefig(overlay, plotsdir("behavior","implant_performance","statespaces_cuememOverlay.png"))

#                               
# ,---.         |,---.|         
# |  _.,---.,---||__. |    ,   .
# |   |,---||   ||    |    |   |
# `---'`---^`---'`    `---'`---|
#                          `---'
Gadfly.push_theme(:dark)
ctab,wtab = tab[tab.wm.==0,:], tab[tab.wm.==1,:]
Gadfly.push_theme(Theme(bar_highlight=nothing))
    coord = Coord.cartesian(xmin=0,xmax=1)
    state  = Gadfly.plot(tab[tab.animal.!="RY9", :], x=:state,    color=:animal, coord, Guide.xlabel("Overall"),       Geom.histogram(position=:identity, density=true));
    cue    = Gadfly.plot(ctab,                       x=:cuestate, color=:animal, coord, Guide.xlabel("Cue-Guided"),    Geom.histogram(position=:identity, density=true));
    mem    = Gadfly.plot(wtab,                       x=:wmstate,  color=:animal, coord, Guide.xlabel("Memory-Guided"), Geom.histogram(position=:identity, density=true))
    stated = Gadfly.plot(tab,                        x=:state,    color=:animal, coord, Guide.xlabel("Overall"),       Theme(alphas=[0.6]), Geom.density);
    cued   = Gadfly.plot(ctab,                       x=:cuestate, color=:animal, coord, Guide.xlabel("Cue-Guided"),    Theme(alphas=[0.6]), Geom.density);
    memd   = Gadfly.plot(wtab,                       x=:wmstate,  color=:animal, coord, Guide.xlabel("Memory-Guided"), Theme(alphas=[0.6]), Geom.density);
    set_default_plot_size(28cm,16cm)
    final = vstack(hstack(state, cue, mem),
           hstack(stated, cued, memd))
final |> SVG(plotsdir("behavior","implant_performance","animal_distribution_summary.svg"))
final |> PDF(plotsdir("behavior","implant_performance","animal_distribution_summary.pdf"))
final |> PNG(plotsdir("behavior","implant_performance","animal_distribution_summary.png"))

ctab,wtab = tab[tab.wm.==0,:], tab[tab.wm.==1,:]
S = style(bar_highlight=nothing, line_width=0cm, highlight_width=0cm)
set_default_plot_size(20cm,16cm)
Gadfly.push_theme(S)

state  = Gadfly.plot(tab[tab.animal.!="RY9", :],x=:id, y=:state,    linestyle=:day, alpha=[0.1], color=:animal, Guide.xlabel("trial"), Guide.ylabel("Performance", orientation=:horizontal), Guide.title("Overall"),       Geom.point);
cue    = Gadfly.plot(ctab,                      x=:id, y=:cuestate, linestyle=:day, alpha=[0.1], color=:animal, Guide.xlabel("trial"), Guide.ylabel("Performance", orientation=:horizontal), Guide.title("Cue-Guided"),    Geom.point);
mem    = Gadfly.plot(wtab,                      x=:id, y=:wmstate,  linestyle=:day, alpha=[0.1], color=:animal, Guide.xlabel("trial"), Guide.ylabel("Performance", orientation=:horizontal), Guide.title("Memory-Guided"), Geom.point);
scat = vstack(state, cue, mem)

tab.id = Vector{Float64}(tab.id)
ctab.id = Vector{Float64}(ctab.id)
wtab.id = Vector{Float64}(wtab.id)
state  = Gadfly.plot(tab[tab.animal.!="RY9", :],x=:id, y=:state,    color=:animal, Guide.xlabel("trial"), Guide.ylabel("Performance", orientation=:vertical), Guide.title("Overall"),       Stat.binmean(n=20), Geom.point, Geom.line)
cue    = Gadfly.plot(ctab,                      x=:id, y=:cuestate, color=:animal, Guide.xlabel("trial"), Guide.ylabel("Performance", orientation=:vertical), Guide.title("Cue-Guided"),    Stat.binmean(n=20), Geom.point, Geom.line);
mem    = Gadfly.plot(wtab,                      x=:id, y=:wmstate,  color=:animal, Guide.xlabel("trial"), Guide.ylabel("Performance", orientation=:vertical), Guide.title("Memory-Guided"), Stat.binmean(n=20), Geom.point, Geom.line);
bmean_sum = vstack(state, cue, mem)
final = hstack(scat, bmean_sum)


final |> SVG(plotsdir("behavior","implant_performance","animal_session_scatterbmean_summary.svg"))
final |> PDF(plotsdir("behavior","implant_performance","animal_session_scatterbmean_summary.pdf"))
final |> PNG(plotsdir("behavior","implant_performance","animal_session_scatterbmean_summary.png"))
bmean_sum |> SVG(plotsdir("behavior","implant_performance","animal_session_bmeansum_summary.svg"))
bmean_sum |> PDF(plotsdir("behavior","implant_performance","animal_session_bmeansum_summary.pdf"))
bmean_sum |> PNG(plotsdir("behavior","implant_performance","animal_session_bmeansum_summary.png"))
scat |> SVG(plotsdir("behavior","implant_performance","animal_session_scatter_summary.svg"))
scat |> PDF(plotsdir("behavior","implant_performance","animal_session_scatter_summary.pdf"))
scat |> PNG(plotsdir("behavior","implant_performance","animal_session_scatter_summary.png"))

#stated = Gadfly.plot(tab,                       x=:id, y=:state,    color=:animal, coord, Guide.ylabel("Overall"),       Geom.line);
#cued   = Gadfly.plot(ctab,                      x=:id, y=:cuestate, color=:animal, coord, Guide.ylabel("Cue-Guided"),    Geom.line);
#memd   = Gadfly.plot(wtab,                      x=:id, y=:wmstate,  color=:animal, coord, Guide.ylabel("Memory-Guided"), Geom.line);
overall = vstack(hstack(state, cue, mem),
       hstack(stated, cued, memd))
overall |> SVG(plotsdir("behavior","implant_performance","overall.svg"))
overall |> PDF(plotsdir("behavior","implant_performance","overall.pdf"))
overall |> PNG(plotsdir("behavior","implant_performance","overall.png"))
