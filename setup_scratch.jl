@time using GoalFetchAnalysis
import DI, DIutils
using DataFrames, Statistics
using Plots
using Infiltrator, ProgressMeter
animal, day = "super_clean", 0
@time SPIKES, BEH, CELLS, RIPPLES = DI.load(animal, day, 
    data_source=["spikes","behavior", "cells", "ripples"])
beh2 = DI.load_behavior("super",0)
SPIKES.cycle = Vector{Union{Float32, Missing}}(missing, size(SPIKES,1)) 
DI.annotate_pyrlayer!(CELLS)

Plot.setparentfolder("isolated")
# subset(combine(groupby(cells, [:animal,:tetrode,:pyrlayer]),
#     :meanrate=>mean,
#     :meanrate=>length=>:n_cells,
#     :propbursts=>mean,
# ), :pyrlayer=>s->s.!=false)
# Time conversion factors
# time_factors = DI.superanimal_timeconversion("super",0)
time_factors  = DI.behavior_0time_before_process("super")
animals, days = unique(CELLS.animal), unique(CELLS.day)
animal, day   = (zip(animals, days)|>collect)[2]
DIutils.pushover("Finished loading ALL ANIMAL data")
oad_from_checkpoint = true

        println("Loading data for $animal, day $day")
        Plot.setappend("_$animal")
        lfp, l_ref = begin 
            lfp = DI.load_lfp(animal, day);
            # lfp = DI.load_lfp(animal, day; vars=["time","tetrode","raw"]);
            # lfp = subset(lfp, :tetrode=>t->t.∈ (unique(l_pyr.tetrode),), 
            #              view=true)
            # DIutils.filtreg.register(groupby(lfp,"tetrode"), 
            #                          groupby(l_pyr,"tetrode"), 
            #                          on="time", transfer=["raw"])
            lfp.time .-= time_factors[animal];
            l_ref =try
                l_ref = DI.load_lfp(animal, day; tet=DI.ca1ref_tetrodes[animal]);
                l_ref.time .-= time_factors[animal];
                l_ref
            catch
                @warn "No reference lfp for $animal"
                l_ref = DataFrame()
            end
            lfp, l_ref
        end
        # Print extrema of spikes.time and lfp.time
        begin
            println("spikes.time: ", extrema(SPIKES.time))
            println("lfp.time: ", extrema(lfp.time))
        end
        
#                                                       
        #  |     ,---.,---.--.--.   .,---.    |    ,---.,---.
        #  |     `---.|---   |  |   ||---'    |    |__. |---'
        #  |         ||      |  |   ||        |    |    |    
        #  `o    `---'`---'  `  `---'`        `---'`    `    
        #                                                    
        # Let's check the consistency of theta across pyr layer tetrodes
        all_tetrodes = pyr_tetrodes = subset(CELLS, :pyrlayer=>s->s.!=false,
            :animal=>a->a.==animal).tetrode |> unique
        all_tetrodes = [all_tetrodes..., DI.default_tetrodes[animal],
                        DI.ca1ref_tetrodes[animal]]
        @time lfp  = subset(lfp, :tetrode=>t->t.∈(all_tetrodes,), view=true);
        GC.gc()

        # Visualize
        # begin
        #     ulfp = unstack(lfp, :time, :tetrode, :raw);
        #     # ulfp = unstack(lfp, :time, :tetrode, :cycle);
        #     sort(Int.(unique(lfp.tetrode)))
         #     Plot.setfolder("lfp")
        #     @time begin
        #         s = 6_000
        #         samp = Matrix(ulfp[s:s+1_000, Not(:time)])
        #         C=cor(Matrix(ulfp[1:100_000, Not(:time)]), dims=1)
        #         xtcks = (1:size(C,1), string.(all_tetrodes))
        #         ytcks = (1:size(C,2), string.(all_tetrodes))
        #         h=heatmap(Matrix(samp), xticks=xtcks, yticks=ytcks, xrotation=90, 
        #             yrotation=0, title="Theta on PYR tetrodes")
        #         p=plot()
        #         least_corr = all_tetrodes[sortperm(mean(C,dims=2)|>vec)][1:3]
        #         most_correlated_tetrode = all_tetrodes[argmax(mean(C,dims=2)|>vec)]
        #         for (col, tet) in zip(eachcol(samp), all_tetrodes)
        #             if tet in pyr_tetrodes && tet != most_correlated_tetrode
        #                 plot!(col, label=tet,  linewidth = (tet in least_corr)  ? 
        #                     3 : 1, linestyle= (tet in least_corr) ? :solid : :dash)
        #             elseif tet == most_correlated_tetrode
        #                 plot!(col, label="M",  linewidth=5, color=:black,
        #                     linestyle=:dash)
        #             else
        #                 plot!(col, label="--", linewidth=1, color=:black,
        #                     linestyle=:dash)
        #             end
        #         end
        #         current()
        #         hc=heatmap(C, xticks=xtcks, yticks=ytcks, xrotation=90, yrotation=0)
        #         layout = @layout [[a b]; c{0.5w}]
        #         plot(h, p, hc, layout=(1,3), size=(1200,400))
        #     end
        #     Plot.save("theta_sync")
        #     ulfp = nothing; 
        #     current()
        # end

        

        l_pyr = subset(lfp, :tetrode=>t->t.∈(pyr_tetrodes,), view=true);
        lfp = nothing; GC.gc()
        # l_def = subset(lfp, :tetrode=>t->t.∈(DI.default_tetrodes[animal],),
        # view=true)
        DIutils.pushover("Ready to ripple - temp!")
        # Get the ripple band
        fs = 1/median(diff(l_pyr.time[1:15_000])) 
        GoalFetchAnalysis.Munge.lfp.bandpass(groupby(l_pyr, :tetrode), 
                                                     150, 250; order=9, fs, 
                                                     newname="ripple")
        l_pyr.raw         = convert(Vector{Union{Missing, Int16}}, 
                                    round.(l_pyr.raw))
        l_pyr.ripple      = convert(Vector{Union{Missing, Float32}}, 
                                    (l_pyr.ripple));
        l_pyr.rippleamp   = convert(Vector{Union{Missing, Float32}},
                                    (l_pyr.rippleamp));
        l_pyr.ripplephase = convert(Vector{Union{Missing, Float32}},
                                    (l_pyr.ripplephase));
        

        DIutils.pushover("Head back to isolated!")

        # ,--.     ,---.     |                                 |              
        # ,--'     `---.,---.|--- .   .,---.    ,---.,   .,---.|    ,---.,---.
        # |            ||---'|    |   ||   |    |    |   ||    |    |---'`---.
        # `--'o    `---'`---'`---'`---'|---'    `---'`---|`---'`---'`---'`---'
        #                            |             `---'                    
        l_pyr[!,:cycle] = Vector{Union{Missing, Int32}}(missing, nrow(l_pyr))
        @time Munge.lfp.annotate_cycles!(groupby(l_pyr, :tetrode), 
                                         method="peak-to-peak")
        # ISSUE: MEMORY explodes on the above line
        l_pyr[!,:cycle] = convert(Vector{Union{Missing, Int32}}, 
                                  l_pyr[!,:cycle]);
        GC.gc()
        # INFO: SAVE
        DI.save_lfp(transform(l_pyr[!,[:time, :tetrode, :cycle, :raw, :phase,
                                    :broadraw, :ripple, :rippleamp, :ripplephase]],
            :time=> t-> t .+ time_factors[animal], renamecols=false), 
            animal, day; handlemissing=:drop, append="pyr")
        DIutils.pushover("Finished annotating cycles for $animal")

        # Examine
        ripples = subset(RIPPLES, :animal=>a->a.==animal, :day=>d->d.==day,
                        view=true)

        using Peaks
        # get second tet
        l = groupby(l_pyr, :tetrode)[2]
        # plot(l.ripple[1:1000], label="ripple",
        #     title="peak count = $(length(Peaks.maxima(l.ripple[1:1000])))")

        i=1000
        ripple_time = ripples.start[i];
        ripple_stop = ripples.stop[i];
        ind  = DIutils.searchsortednearest(l.time, ripple_time);
        ind1 = DIutils.searchsortednearest(l.time, ripple_stop);
        plot(l.ripple[ind:ind1], label="ripple",
            title="peak count = $(length(Peaks.maxima(l.ripple[ind:ind+1000])))")
        plot!(l.ripplephase[ind:ind1]*maximum(l.ripple[ind:ind1]), label="ripple phase")
        plot!(l.rippleamp[ind:ind1], label="ripple amp")
        plot!(DIutils.norm_extrema(l.broadraw[ind:ind1]).*maximum(l.rippleamp[ind:ind1]), label="raw")

        # BUG: ERROR IN THIS SECTION
        # Get cycle table
        cycles = Munge.lfp.get_cycle_table(groupby(l_pyr, :tetrode));
        cycles.cycle = convert(Vector{Union{Missing, Int32}}, cycles.cycle);
        for col in [:start, :stop, :amp_mean, :δ]
            cycles[!,col] = convert(Vector{Union{Missing, Float32}}, 
                cycles[!,col])
        end
        # INFO: SAVE
        # Save cycles
        DI.save_cycles(begin
            transform(cycles, 
                :start=>s->s .+ time_factors[animal],
                :stop=>s->s  .+ time_factors[animal],
                renamecols=false
            )
        end, animal, day, "pyr")
        DIutils.pushover("Finished saving cycles for $animal")

        # ---------------- *
        # Visualize cycles |
        # ---------------- *
        using GoalFetchAnalysis.Plot.lfplot
        begin
            folder, par = Plot.folder_args, Plot.parent_folder
            Plot.setparentfolder("theta");
            Plot.setfolder("cycle_cutting");
            P = []
            for g in groupby(l_pyr, :tetrode)
                push!(P,lfplot.cycleplot(g))
            end
            plot(P[1:10]..., layout=(10,1), size=(1200, 400*10))
            Plot.save("cycleplot_of_pyr_tetrodes")
            # Restore settings
            Plot.setfolder(folder...) 
            Plot.setparentfolder(par...)
        end

    # ISSUE: ONE USUALLY HAS TO RELOAD FROM CHECKPOINTED FILES HERE UNLESS
    # YOU HAVE ENORMOUS RAM partially because of the lfp
    # CONSIDER ---> reducing col sizes
    #
    #                                                                 
    # ,--.     ,---.     |                            o|              
    #   -|     `---.,---.|--- .   .,---.    ,---.,---..|__/ ,---.,---.
    #    |         ||---'|    |   ||   |    `---.|   |||  \ |---'`---.
    # `--'o    `---'`---'`---'`---'|---'    `---'|---'``   ``---'`---'
    #                              |             |                    
    ripple_props = [:ripple, :ripple_time, :ripple_phase, :ripple_phase_band,
                    :ripple_amp_band];
    # Ripple phase
    spikes  = subset(SPIKES,  :animal=>a->a.==animal, :day=>d->d.== day,
                    view=true);
    ripples = subset(RIPPLES, :animal=>a->a.==animal, :day=>d->d.== day,
                    view=true);
    Munge.spiking.event_spikestats!(spikes, ripples; eventname="ripple");
    @assert(length(unique(spikes.ripple_phase)) .> 4, 
    "Ripple phase not calculated correctly")

    # ISSUE:
    # Let's figure out the number of tetrodes that match between
    # spikes and lfp
    sp = subset(spikes, 
    :tetrode=>t->t .∈ (unique(l_pyr.tetrode),), view=true)
    println("Theoretical number of pyramidal layer cells:",
    combine(groupby(sp, :tetrode),
    :unit => (x->length(unique(x))), renamecols=false)
    )

    # Add ripple BAND phase (rather than phase within the ripple)
    spikes.ripple_phase_band = Vector{Union{Missing, Float32}}(missing, 
        size(spikes,1));
    spikes.ripple_amp_band   = Vector{Union{Missing, Float32}}(missing, 
        size(spikes,1));
    gl_pyr   = groupby(l_pyr,  :tetrode);
    g_spikes = groupby(spikes, :tetrode);
    K = intersect(keys(gl_pyr)  |>collect .|> NamedTuple,
        keys(g_spikes)|>collect           .|> NamedTuple);
    println("Number of matching keys:", length(K))

    @assert all(ripple_props .∈ (propertynames(spikes),))
    using ProgressMeter
    prog = Progress(length(K), 1, "Adding ripple band phase");
    k = (;tetrode = 55)
    for k in K
        gs, gl = g_spikes[k], gl_pyr[k]
        time = gs.time
        ind = DIutils.searchsortednearest.([gl.time], time)
        gs.ripple_phase_band = gl.ripplephase[ind];
        gs.ripple_amp_band   = gl.rippleamp[ind];
        next!(prog)
    end

    # ISSUE: Let's figure out what fraction of our pyramidal layer
    #        cells
    q=combine(groupby(sp, :tetrode),
    :unit => (x->length(unique(x))), 
    [:unit, :ripple_phase_band] => ((x,y)->( goodinds = .!ismissing.(y); unique(x[goodinds]) |> length)) => :n_with_ripple_phase,
        renamecols=false,
    )
    println("Theoretical number of pyramidal layer cells:", q)
    

    # Theta phase
    # sort!(spikes, [:unit, :time])
    # cells = subset(CELLS, :animal=>a->a.==animal, :day=>d->d.==day,
    #                 view=true)
    # DIutils.filtreg.register(cells, spikes, on="unit", transfer=["tetrode"])
    Munge.spiking.event_spikestats!(groupby(spikes,:tetrode), 
                                    groupby(cycles, :tetrode); eventname="theta")
    theta_props = [:theta, :theta_time, :theta_phase]
    @assert all(theta_props .∈ (propertynames(spikes),))

    # PREPARE FOR ISOLATED SPIKE ANALYSIS
    GC.gc()
    # Cycles filtered by ripples and movement
    spikes=GoalFetchAnalysis.Munge.spiking.prepiso(
        spikes, l_pyr; 
        cycle=:cycle,
        cells=CELLS, include_samples=false, matchtetrode=true, 
        BEH, RIPPLES)
    # Cycles unfiltered
    spikes.unfilt_cycle = Vector{Union{Missing,Bool}}(missing, nrow(spikes))
    l_pyr[!,:unfilt_cycle] = l_pyr[!, :cycle]
    spikes=GoalFetchAnalysis.Munge.spiking.prepiso(
        spikes, l_pyr; 
        cycle=:unfilt_cycle, cells=CELLS, include_samples=false,
        matchtetrode=true)
    # Sanity check : should be more missing cycles in the filtered vs 
    #   unfiltered
    #   
    cycle_props = [:animal, :day, :time, :tetrode, :unit, :cycle, 
                   :unfilt_cycle]
    props = intersect(union(cycle_props, ripple_props), propertynames(spikes))
    DI.save_spikes_taginfo(begin
        transform(spikes[!,props], :time=>t->t.+time_factors[animal], renamecols=false)
    end, animal, day, "pyr_cycles_unfilt")
    println("Missing cycles in filtered: ",   sum(ismissing.(spikes.cycle)))
    println("Missing cycles in unfiltered: ", sum(ismissing.(spikes.unfilt_cycle)))

    # Sanity check : should be more missing cycles in the filtered vs
    @assert(sum(ismissing.(spikes.unfilt_cycle)) > 
            sum(ismissing.(spikes.cycle)),
            "More missing cycles in filtered than unfiltered")

    # Setup isolated spiking fields
    if :isolated ∉ propertynames(spikes)
        spikes.isolated  = Vector{Union{Bool, Missing}}(missing, nrow(spikes))
    end
    sp = nothing
    using ProfileView
    # INFO: APPARENTLY runs faster if function in Main
    # @profile begin
        sp = Munge.spiking.isolated!(spikes)
    # end
    # ProvileView.view()
    # ISSUE: invesetigate long units like 134: cell has a meanrate=3.5hz)
    iso_props = [:isolated, :nearestcyc, :meancyc]
    props = intersect(union(cycle_props, ripple_props, theta_props, iso_props),
                      propertynames(spikes))
    DI.save_spikes_taginfo(begin
        transform(spikes[!,props], 
            :time=>t->t.+time_factors[animal],
        renamecols=false)
    end, animal, day, "pyr_cycles_isolated")

