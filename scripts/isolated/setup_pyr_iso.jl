@time using GoalFetchAnalysis
import DI, DIutils
using DataFrames, Statistics
using Plots
using Infiltrator
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
animal, day   = zip(animals, days)|>first
DIutils.pushover("Finished loading ALL ANIMAL data")
animal = animals[1]
load_from_checkpoint = true

for animal in animals

    # CHECKPOINT LOAD
    if load_from_checkpoint && isfile(DI.cyclepath(animal,day,"pyr"))

        l_pyr = DI.load_lfp(animal, day; append="pyr")
        l_pyr = transform(l_pyr, :time=> t-> t .- time_factors[animal], renamecols=false)
        # transform!(l_pyr, :time=> t-> t .+ time_factors[animal], renamecols=false)
        cycles = DI.load_cycles(animal, day, "pyr")
        cycles.start = cycles.start_function
        cycles.stop  = cycles.stop_function
        cycles = transform(cycles,
            :start=> t-> t .- time_factors[animal],
            :stop=> t-> t .- time_factors[animal],
            renamecols=false
        )
        cycles = cycles[!,Not([:start_function,:stop_function])]
        spikes = transform(DI.load_spikes(animal, day, "pyr_cycles_isolated"),
            :time=> t-> t .- time_factors[animal],
            renamecols=false
        )
        # spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day)
        # Print out extrema(time) for each of these just to be sure
        
        begin
            println("spikes.time: ", extrema(spikes.time))
            println("lfp.time: ", extrema(l_pyr.time))
            println("cycles.time: ", extrema(cycles.start))
        end
        DIutils.pushover("Loaded checkpoint for $animal")

    # CHECKPOINT GENERATE
    else

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
        begin
            ulfp = unstack(lfp, :time, :tetrode, :raw);
            # ulfp = unstack(lfp, :time, :tetrode, :cycle);
            sort(Int.(unique(lfp.tetrode)))
            Plot.setfolder("lfp")
            @time begin
                s = 6_000
                samp = Matrix(ulfp[s:s+1_000, Not(:time)])
                C=cor(Matrix(ulfp[1:100_000, Not(:time)]), dims=1)
                xtcks = (1:size(C,1), string.(all_tetrodes))
                ytcks = (1:size(C,2), string.(all_tetrodes))
                h=heatmap(Matrix(samp), xticks=xtcks, yticks=ytcks, xrotation=90, 
                    yrotation=0, title="Theta on PYR tetrodes")
                p=plot()
                least_corr = all_tetrodes[sortperm(mean(C,dims=2)|>vec)][1:3]
                most_correlated_tetrode = all_tetrodes[argmax(mean(C,dims=2)|>vec)]
                for (col, tet) in zip(eachcol(samp), all_tetrodes)
                    if tet in pyr_tetrodes && tet != most_correlated_tetrode
                        plot!(col, label=tet,  linewidth = (tet in least_corr)  ? 
                            3 : 1, linestyle= (tet in least_corr) ? :solid : :dash)
                    elseif tet == most_correlated_tetrode
                        plot!(col, label="M",  linewidth=5, color=:black,
                            linestyle=:dash)
                    else
                        plot!(col, label="--", linewidth=1, color=:black,
                            linestyle=:dash)
                    end
                end
                current()
                hc=heatmap(C, xticks=xtcks, yticks=ytcks, xrotation=90, yrotation=0)
                layout = @layout [[a b]; c{0.5w}]
                plot(h, p, hc, layout=(1,3), size=(1200,400))
            end
            Plot.save("theta_sync")
            ulfp = nothing; 
            current()
        end

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

        l_pyr.raw = convert(Vector{Union{Missing, Int16}}, round.(l_pyr.raw))
        l_pyr.ripple = convert(Vector{Union{Missing, Float32}}, (l_pyr.ripple))
        l_pyr.rippleamp = convert(Vector{Union{Missing, Float32}}, 
                                  (l_pyr.rippleamp))
        l_pyr.ripplephase = convert(Vector{Union{Missing, Float32}}, 
                                    (l_pyr.ripplephase))
        

        # ,--.     ,---.     |                                 |              
        # ,--'     `---.,---.|--- .   .,---.    ,---.,   .,---.|    ,---.,---.
        # |            ||---'|    |   ||   |    |    |   ||    |    |---'`---.
        # `--'o    `---'`---'`---'`---'|---'    `---'`---|`---'`---'`---'`---'
        #                            |             `---'                    
        @time Munge.lfp.annotate_cycles!(groupby(l_pyr, :tetrode), 
                                         method="peak-to-peak")
        # ISSUE: MEMORY explodes on the above line
        l_pyr[!,:cycle] = convert(Vector{Union{Missing, Int32}}, 
                                  l_pyr[!,cycle]);
        GC.gc()
        DI.save_lfp(transform(l_pyr[!,[:time, :tetrode, :cycle, :raw, :phase,
                                    :broadraw, :ripple, :rippleamp, :ripplephase]],
            :time=> t-> t .+ time_factors[animal], renamecols=false), 
            animal, day; handlemissing=:drop, append="pyr")
        DIutils.pushover("Finished annotating cycles for $animal")

        # Get cycle table
        cycles = Munge.lfp.get_cycle_table(groupby(l_pyr, :tetrode))
        cycles.cycle = convert(Vector{Union{Missing, Int32}}, cycles.cycle)
        for col in [:start, :stop, :amp_mean, :δ]
            cycles[!,col] = convert(Vector{Union{Missing, Float32}}, 
                cycles[!,col])
        end
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
        begin
            folder, parent = Plot.folder_args, Plot.parent_folder
            Plot.setparentfolder("theta");
            Plot.setfolder("cycle_cutting")
            using GoalFetchAnalysis.Plot.lfplot
            P = []
            for g in groupby(l_pyr, :tetrode)
                push!(P,lfplot.cycleplot(g))
            end
            plot(P[1:10]..., layout=(10,1), size=(1200, 400*10))
            Plot.save("cycleplot_of_pyr_tetrodes")
            # Restore settings
            Plot.setfolder(folder...) Plot.setparentfolder(parent...)
        end


    end

    #                                                                 
    # ,--.     ,---.     |                            o|              
    #   -|     `---.,---.|--- .   .,---.    ,---.,---..|__/ ,---.,---.
    #    |         ||---'|    |   ||   |    `---.|   |||  \ |---'`---.
    # `--'o    `---'`---'`---'`---'|---'    `---'|---'``   ``---'`---'
    #                              |             |                    
    
    
    # Ripple phase
    spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day,
                    view=true)
    ripples = subset(RIPPLES, :animal=>a->a.==animal, :day=>d->d.==day,
                    view=true)
    DIutils.pushover("Ready to check RIPs")
    Munge.spiking.event_spikestats!(spikes, ripples; eventname="ripple")
    ripple_props = [:ripple, :ripple_time, :ripple_phase]
    @assert all(ripple_props .∈ (propertynames(spikes),))
    @assert(length(unique(spikes.ripple_phase)) .> 4, "Ripple phase not"* 
        "calculated correctly")

    # Add ripple BAND phase (rather than phase within the ripple)
    spikes.ripple_phase_band = Vector{Union{Missing, Float32}}(missing, 
        size(spikes,1))
    gl_pyr = groupby(l_pyr, :tetrode)
    prog = Progress(length(spikes), 1, "Adding ripple band phase")
    Threads.@threads for spike in eachrow(spikes)
        time =spike.time
        tetrode = spike.tetrode
        l = gl_pyr[(;tetrode)]
        ind = DIutils.searchsortednearest(l.time, time)
        spike.ripple_phase_band = l.ripplephase[ind]
        next!(prog)
    end

    # Theta phase
    # sort!(spikes, [:unit, :time])
    # cells = subset(CELLS, :animal=>a->a.==animal, :day=>d->d.==day,
    #                 view=true)
    # DIutils.filtreg.register(cells, spikes, on="unit", transfer=["tetrode"])
    Munge.spiking.event_spikestats!(groupby(spikes,:tetrode), 
                                    groupby(cycles, :tetrode); eventname="theta")
    theta_props = [:theta, :theta_time, :theta_phase]
    @assert all(theta_props .∈ (propertynames(spikes),))

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
    cycle_props = [:animal, :day, :time, :tetrode, :unit, :cycle, :unfilt_cycle]
    props = intersect(union(cycle_props, ripple_props), propertynames(spikes))
    DI.save_spikes_taginfo(begin
        transform(spikes[!,props], 
        :time=>t->t.+time_factors[animal],
    renamecols=false)
    end, animal, day, "pyr_cycles_unfilt")
    println("Missing cycles in filtered: ",   sum(ismissing.(spikes.cycle)))
    println("Missing cycles in unfiltered: ", sum(ismissing.(spikes.unfilt_cycle)))
    @assert sum(ismissing.(spikes.cycle)) >   sum(ismissing.(spikes.unfilt_cycle))

    # Setup isolated spiking fields
    spikes.isolated = Vector{Union{Bool, Missing}}(missing, nrow(spikes))
    sp = nothing
    using ProfileView
    @profile begin
    sp = Munge.spiking.isolated!(spikes)
    end
    ProvileView.view()
    # ISSUE: invesetigate long units like 134: cell has a meanrate=3.5hz)
    iso_props = [:isolated, :nearestcyc, :meancyc]
    props = intersect(union(cycle_props, ripple_props, theta_props, iso_props),
                      propertynames(spikes))
    DI.save_spikes_taginfo(begin
        transform(spikes[!,props], 
            :time=>t->t.+time_factors[animal],
        renamecols=false)
    end, animal, day, "pyr_cycles_isolated")


end

