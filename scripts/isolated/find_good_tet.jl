using DrWatson
quickactivate(expanduser("~/Projects/goal-code/"));
using GoalFetchAnalysis
using DIutils.namedtup
Munge = GoalFetchAnalysis.Munge
Filt = GoalFetchAnalysis.Filt
Plot = GoalFetchAnalysis.Plot
import GoalFetchAnalysis.Munge: spiking
import DI, DIutils

using Infiltrator, SoftGlobalScope
using DataStructures: OrderedDict
using DimensionalData
import DimensionalData: Between
using ProgressMeter
using DataFrames, DataFramesMeta
using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
using Plots
using LazySets
import FileIO
using LaTeXStrings
# using AlgebraOfGraphics
# set_aog_theme!()
# AoG = AlgebraOfGraphics

datasets = DI.animaldays()
filt = Filt.get_filters()
datacut = :all
results = OrderedDict()
(animal, day) = first(datasets)
spikes = nothing



function add_suptitle(P, suptitle)
    # Create a subplot with the title but no axes
    suptitle_subplot = Plots.plot(legend=false, grid=false)
    Plots.plot!(title=suptitle, framestyle=:none)

    if P isa Vector
        P = Plots.plot(P...)
    end

    # Combine the suptitle plot with the original plot
    P_new = Plots.plot(suptitle_subplot, P, 
        layout=@layout([a{0.025h}; b{0.975h}]), 
        margin=0Plots.mm)

    return P_new
end


        function showtimes()
        println("animal=$animal, day=$day, tet=$tet",
            "\nlfp: ", extrema(lfp.time), 
            "\nspikes: ", extrema(spikes.time),
            "\nbeh: ", extrema(beh.time),
            "\nripples: ", extrema(ripples.start), " ", extrema(ripples.stop)
        )
        end


@showprogress "animal" for (animal,day) in datasets

    global spikes, lfp

    Plot.setappend((;animal,day))

    # ===================
    # ACQUIRE DATA
    # ===================
    # Acquire data
    @time spikes, beh, cells, ripples = DI.load(animal, day, 
        data_source=["spikes","behavior", "cells", "ripples"])
    GC.gc()
    # sort!(spikes, :time)
    # sort!(beh, :time)
    # @assert typeof(spikes) <: DataFrame
    # @assert typeof(beh) <: DataFrame
    # beh, spikes = DIutils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], 
    # filters=filt[datacut], filter_skipmissingcols=true)
    # allspikes = copy(spikes)
    # beh2 = DI.load_behavior(animal,day)
    # Munge.nonlocal.setunfilteredbeh(beh2)

    DI.annotate_interneuron!(cells)
    ca1 = sort(unique(@subset(cells, :area .== "CA1").tetrode))
    pyr_units = @subset(cells, :area .== "CA1", :interneuron .== false).unit
    GC.gc()

    spikes = subset(spikes, :area=>a->a .== "CA1", 
        :unit=>u->u .∈ (pyr_units,), view=true)


    iso, ncells, npyr = OrderedDict(), OrderedDict(), OrderedDict(),
                                OrderedDict()
    ripmeans, rips = OrderedDict(), OrderedDict()
    tet = first(ca1)

    @showprogress "gathering" for tet in ca1#= |>x->Iterators.take(x,2) =#

        global spikes, lfp
        println("Tetrode $tet")
        println("Exists spikes var :", isdefined(Main,:spikes))
        lfp = DI.load_lfp(animal, day; tet=tet,
            vars=["time","raw","amp","phase", "broadraw"])
        lfp.time = lfp.time .- DI.min_time_records[end]


        if :phase in propertynames(spikes)
            spikes = DataFrame(spikes[!, Not(:phase)], copycols=false)
        end
        if :cycle in propertynames(spikes)
            spikes = DataFrame(spikes[!, Not(:cycle)], copycols=false)
        end
        # try
        Munge.lfp.annotate_cycles(lfp, method="peak-to-peak") # TODO potential bug, 1st time runs, cuts trough-to-trough, second peak-to-peak
        Munge.lfp.handle_phase_asymmetry!(lfp, phase_col_out="corrected_phase")
        # catch
        #     push!(iso, tet=>nothing)
        #     continue
        # end
        @assert length(unique(lfp.cycle)) > 1
        #sp = @subset(spikes, :tetrode .== 6);
        # Does the cell ceom from this tetrode?
        spikes[!,:match_tet] = spikes.tetrode .== tet
        # Make phase less ram expenseive
        lfp.phase = convert(Array{Union{Missing,Float32}}, lfp.phase)
        # Add phases to spikes
        println("Phase locking for tetrode $tet")

        print("...sorting")
        sort!(lfp, :time)
        sort!(spikes, :time)
        print("...registering")
        lfp, spikes = DIutils.filtreg.register(lfp, spikes, on="time", 
                    transfer=["phase", "amp","cycle", "corrected_phase"])
        try
            @assert !all(ismissing.(spikes.phase))
            @assert !all(ismissing.(spikes.cycle))
        catch
            @infiltrate
        end
        # Cells
        # 1. how many cells?
        println("Cell count for tetrode $tet")
        push!(ncells, tet=>length(unique(@subset(cells, :tetrode .== tet).unit)))
        push!(npyr,   tet=> length(unique(@subset(cells, :tetrode .== tet, 
            :interneuron .== false).unit)))

        # Ripples
        # Let's grab each ripples for this tetrode
        Munge.lfp.bandpass(lfp, 150, 250, order=6, newname="ripple")
        Munge.lfp.bandpass(lfp,  0.5,  2, order=3, newname="sharpwave")
        push!(ripmeans, tet=>describe(lfp[!, [:rippleamp, :sharpwaveamp, :ripplephase, :sharpwavephase, :sharpwave, :ripple]]))

        push!(rips, tet=>DIutils.in_range(lfp, ripples))

        # ===================
        # ISOLATED SPIKING
        # ===================
        #Munge.spiking.isolated(sp, lfp, include_samples=true)
        println("Isolated spiking for tetrode $tet")
        Munge.spiking.prepiso(
            spikes, lfp; 
            cycle=:cycle,
            cells=cells, include_samples=false, refreshcyc=false,
            beh, ripples)
        if :isolated ∉ propertynames(spikes)
            spikes.isolated  = Vector{Union{Bool, Missing}}(missing, nrow(spikes))
        end
        # spikes = Munge.spiking.isolated(spikes)
        push!(iso, tet=>spikes[!,[:time,:phase,:amp,:corrected_phase,:match_tet,#=:isolated=#]])
    end
    push!(results, animal=>(iso=iso, ncells=ncells, npyr=npyr, rips=rips, ripmeans=ripmeans))

end

#: SECTION: READY TO PLOT

"""
Function that takes the describe output of a ripple and sharpwave and expands it into a dict where each field is the :variable column item and its value is the corresponding DataFrameRow
"""
function expand_ripmeans(ripmeans)
    ripmeans = ripmeans[!, [:variable, :mean, :min, :max, :nmissing, :eltype, :median]]
    ripmeans = ripmeans[completecases(ripmeans), :]
    variables = unique(ripmeans.variable)
    results = Dict()
    for var in variables
        results[var] = ripmeans[findfirst(ripmeans.variable.==var), :]
    end
    return results
end

file = datadir("checkpoint__find_good_tet.jl.tmp.jld2")
FileIO.save(file, Dict("results"=>results); compress=true)
results = FileIO.load(file)["results"]

DF = []
for (animal,V) in results
    tets = keys(V.iso)
    for tet in tets
        println("Animal $animal, tetrode $tet")
        df = DataFrame(V.iso[tet], copycols=false)
        df.tetrode = fill(tet, nrow(df))
        df.animal = fill(animal, nrow(df))
        df.ncells = fill(V.ncells[tet], nrow(df))
        df.npyr = fill(V.npyr[tet], nrow(df))
        push!(DF, df)
    end
end
DF = vcat(DF...)


#: SECTION: PLOTS w/o AoG
using ElectronDisplay
Plot.setparentfolder("isolated","find_goot_tet.jl")
Plot.theme(:dracula)

animal, day = datasets[2]
DF = subset(DF, :animal=>a->a .== animal, view=true)
Plot.setappend("$(animal)_$(day)")
ripmeans = Dict(k=>expand_ripmeans(v) 
    for (k,v) in results[animal].ripmeans)

gettitle(df, tet) = 
    replace("$(tet), ncells=$(df.ncells[1]), npyr=$(df.npyr[1])"* 
    "\nspw median = $(ripmeans[tet][:sharpwave].median)", 
    "ncells" => L"$n_{cells}$",
    "npyr" => L"$n_{pyr}$", )

# SUBSECTION: Phase locking for all spikes to each tetrode
P = []
ca1=DF.tetrode|>unique
tet=first(ca1)
for tet in ca1
    df = subset(DF, :tetrode=>t-> t.== tet, view=true)
    pl = Plots.histogram(df.phase, bins=40;
        fill=0,
        normalize=:probability,label="",xticks=[],yticks=[], alpha=0.4,
        title=gettitle(df, tet)
    )
    # get max and min of histogram bins
    Plots.ylims!(extrema(pl.series_list[1][:y]))
    push!(P, pl)
end
PP=add_suptitle(
Plots.plot(P...; titlefontsize=6, size=(800,800)),
    "\n\n ALL"
)
Plot.save("all")
current()

# SUBSECTION: pl only for different tetrodes
P,phases,amps = [],[],[]
tet=first(ca1)
for tet in DF.tetrode|>unique
    println("Tetrode $tet")
    df = subset(DF, :tetrode=>t-> t.== tet, 
        :match_tet=>t-> t.== false,
        view=true)
    if !isempty(df)
        pl = Plots.histogram(df.phase, bins=40;
            fill=0,
            normalize=:probability,label="",xticks=[],yticks=[], alpha=0.4,
            title=gettitle(df, tet)
        )
        # get max and min of histogram bins
        Plots.ylims!(extrema(pl.series_list[1][:y]))
        push!(P, pl)
        push!(phases, df.phase)
        push!(amps, df.amp)
    end
end
PP=add_suptitle(
Plots.plot(P...; titlefontsize=6, size=(800,800)),
    "\n\n\n Non-match"
)
Plot.save("nonmatch")
current()

# SUBSECTION: pl only for same tetrode
P,phases,amps = [],[],[]
tet=first(ca1)
for tet in DF.tetrode|>unique
    println("Tetrode $tet")
    df = subset(DF, :tetrode=>t-> t.== tet, 
        :match_tet=>t-> t.== true,
        view=true)
    if !isempty(df)
        pl = Plots.histogram(df.phase, bins=40;
            fill=0,
            normalize=:probability,label="",xticks=[],yticks=[], alpha=0.4,
            title=gettitle(df, tet)
        )
        # get max and min of histogram bins
        Plots.ylims!(extrema(pl.series_list[1][:y]))
        push!(P, pl)
        push!(phases, df.phase)
        push!(amps, df.amp)
    end
end
PP=add_suptitle(
Plots.plot(P...; titlefontsize=6, size=(800,800)),
    "\n\n\n\n\n Match"
)
Plot.save("match")
current()


# SECTION : AMPLITUDE FILTRATION
Plot.setfolder("amp-filtered_theta_phase_locking")
histogram(DF.amp, group=DF.tetrode)
Plot.save("amp-hist")
current()

P = []
for (i,tet) in ca1|>enumerate
    df = subset(DF, :tetrode=>t-> t.== tet, view=true)
    pl = Plots.histogram(df.amp, bins=40, c=:skyblue, fill=0)
    q = quantile(df.amp, [0.1, 0.95, 0.99, 0.999])
    Plots.vline!(q[1:3], color=:white)
    xlims!(0, q[4])
    push!(P, pl)
end
Plots.plot(P...; titlefontsize=6, size=(1200,600), label="")
Plot.save("amp-hist_split_with_quantiles")


# SUBSECTION: pl only for same tetrode
P,phases,amps = [],[],[]
tet=first(ca1)
for tet in DF.tetrode|>unique
    println("Tetrode $tet")
    df = subset(DF, :tetrode=>t-> t.== tet, 
        :match_tet=>t-> t.== true,
        view=true)
    if isempty(df)
        continue
    end
    q = quantile(df.amp, [0.1,0.99])
    df = subset(df, :amp=>a-> a .> q[1] .&& a .<= q[2], view=true)
    if !isempty(df)
        pl = Plots.histogram(df.phase, bins=40;
            fill=0,
            normalize=:probability,label="",xticks=[],yticks=[], alpha=0.4,
            title=gettitle(df, tet)
        )
        # get max and min of histogram bins
        Plots.ylims!(extrema(pl.series_list[1][:y]))
        push!(P, pl)
    end
end
PP=add_suptitle(
Plots.plot(P...; titlefontsize=6, size=(800,800)),
    "\n\n\n\n\n Match"
)
Plot.save("match_amp_filtered")
current()

# SUBSECTION: pl only for same tetrode
P,phases,amps = [],[],[]
tet=first(ca1)
for tet in DF.tetrode|>unique
    println("Tetrode $tet")
    df = subset(DF, :tetrode=>t-> t.== tet, 
        :match_tet=>t-> t.== false,
        view=true)
    if isempty(df)
        continue
    end
    q = quantile(df.amp, [0.1,0.99])
    df = subset(df, :amp=>a-> a .> q[1] .&& a .<= q[2], view=true)
    if !isempty(df)
        pl = Plots.histogram(df.phase, bins=40;
            fill=0,
            normalize=:probability,label="",xticks=[],yticks=[], alpha=0.4,
            title=gettitle(df, tet)
        )
        # get max and min of histogram bins
        Plots.ylims!(extrema(pl.series_list[1][:y]))
        push!(P, pl)
    end
end
PP=add_suptitle(
Plots.plot(P...; titlefontsize=6, size=(800,800)),
    "\n\n NonMatch"
)
Plot.save("NONmatch_amp_filtered")
current()

# # SECTION: PLOTS
#
# using AlgebraOfGraphics
# AoG = AlgebraOfGraphics
# using DataFrames
# using GLMakie
#
# # Assuming DF is your dataframe
# # Filter dataframe for the first animal
# first_animal = unique(DF.animal)[1]
# DF_filtered = DF[DF.animal .== first_animal, :]
#
# # Limit the data to the first 20 tetrodes
# tetrodes = unique(DF_filtered.tetrode)
# if length(tetrodes) > 20
#     tetrodes = tetrodes[1:20]
# end
# DF_filtered = DF_filtered[DF_filtered.tetrode .∈ (tetrodes,), :]
# DF_filt = Dict{Symbol}{Vector}()
# for (key, value) in zip(propertynames(DF_filtered), eachcol(DF_filtered))
#     DF_filt[key] = disallowmissing(value)
# end
# DF_filt = NamedTuple(DF_filt)
#
# # Specify the graphical mapping
# specs = data(DF_filt) * mapping(:phase; layout=:tetrode) * AoG.histogram(bins=range(0,2pi,length=15))
# # Draw the plot
# draw(specs)
#
#
#
#
# # Add the title with the ncells and npyr count
# plt = plt * mapping(layout_y = :ncells, title = :npyr)
#
# # Configure the layout
# plt = plt |> layout_wrap(10)
#
#
#
#
# df = (x=randn(5000), y=randn(5000), z=rand(["a", "b", "c"], 5000))
# specs = data(df) * mapping(:x; layout=:z) * AoG.histogram(bins=range(-2, 2, length=15))
# draw(specs)
#
# # gdf = groupby(DF, :tetrode)
# # pl = [(@df DataFrame(dropmissing(g, :isolated)) histogram(:phase,
# #     group=:isolated, bins=40;
# #     normalize=:probability,label="",xticks=[],yticks=[], alpha=0.4,
# #     title="$(g.tetrode[1])")) for g in gdf] plot(pl...)
# #
# # Plot.setparentfolder("isolated","find_goot_tet.jl")
# # Plot.save("tetrodePhaseLocking")
