@time using GoalFetchAnalysis
using DataFrames, Statistics
using Plots
animal, day = "super_clean", 0

@time SPIKES, BEH, CELLS, RIPPLES = DI.load(animal, day, 
    data_source=["spikes","behavior", "cells", "ripples"])
beh2 = DI.load_behavior("super",0)
SPIKES.cycle = Vector{Union{Float32, Missing}}(missing, size(SPIKES,1)) 
DI.annotate_pyrlayer!(CELLS)
Plot.setparentfolder("isolated")
# 
# subset(combine(groupby(cells, [:animal,:tetrode,:pyrlayer]),
#     :meanrate=>mean,
#     :meanrate=>length=>:n_cells,
#     :propbursts=>mean,
# ), :pyrlayer=>s->s.!=false)


# Time conversion factors
# time_factors = DI.superanimal_timeconversion("super",0)
time_factors = DI.behavior_0time_before_process("super")
animals, days = unique(CELLS.animal), unique(CELLS.day)
animal, day = zip(animals, days)|>first
DIutils.pushover("Finished loading non-lfp")
animal = animals[1]

load_from_checkpoint = true

for animal in animals

    Plot.setappend("_$animal")
    lfp, l_ref = begin 
        lfp = DI.load_lfp(animal, day);
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

    # CHECKPOINT LOAD
    if load_from_checkpoint && isfile(DI.cyclepath(animal,day,"pyr"))

        l_pyr = DI.load_lfp(animal, day; append="pyr")
        l_pyr = 
            transform(l_pyr, :time=> t-> t .- time_factors[animal],
                      renamecols=false)
        cycles = DI.load_cycles(animal, day, "pyr")
        cycles.start = cycles.start_function
        cycles.stop  = cycles.stop_function
        cycles = transform(cycles,
            :start=> t-> t .- time_factors[animal],
            :stop=> t-> t .- time_factors[animal],
            renamecols=false
        )
        cycles = cycles[!,Not([:start_function,:stop_function])]
        spikes = transform(DI.load_spikes(animal, day, "pyr_cycles"),
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

    # CHECKPOINT GENERATE
    else

        # Let's check the consistency of theta across pyr layer tetrodes
        all_tetrodes = pyr_tetrodes = subset(CELLS, :pyrlayer=>s->s.!=false,
            :animal=>a->a.==animal).tetrode |> unique
        all_tetrodes = [all_tetrodes..., DI.default_tetrodes[animal],
                        DI.ca1ref_tetrodes[animal]]
        @time lfp  = subset(lfp, :tetrode=>t->t.∈(all_tetrodes,), view=true);
        ulfp = unstack(lfp, :time, :tetrode, :raw);
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
                    plot!(col, label=tet,  linewidth = (tet in least_corr)  ? 3 : 1,
                        linestyle= (tet in least_corr) ? :solid : :dash)
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

        l_pyr = subset(lfp, :tetrode=>t->t.∈(pyr_tetrodes,), view=true);
        lfp = nothing; GC.gc()
        # l_def = subset(lfp, :tetrode=>t->t.∈(DI.default_tetrodes[animal],),
        # view=true)
        
        # Setup cycles
        @time Munge.lfp.annotate_cycles!(groupby(l_pyr, :tetrode), 
                                         method="peak-to-peak")
        # ISSUE: MEMORY explodes on the above line
        l_pyr[!,:cycle] = convert(Vector{Union{Missing, Int32}}, 
                                  l_pyr[!,cycle]);
        GC.gc()
        DI.save_lfp(transform(l_pyr[!,[:time, :tetrode, :cycle]],
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

    end

    # α Setup spikes
    spikes = subset(SPIKES, :animal=>a->a.==animal, :day=>d->d.==day,
                    view=true)
    GC.gc()
    # Cycles filtered by ripples and movement
    spikes=GoalFetchAnalysis.Munge.spiking.prepiso(
        spikes, l_pyr; 
        cycle=:cycle,
        cells=CELLS, include_samples=false, matchtetrode=true, 
        BEH, RIPPLES)
    DI.save_spikes_taginfo(begin
        transform(spikes[!,[:animal, :day, :time, :unit, :cycle]], :time=>t->t.+time_factors[animal]
        renamecols=false)
    end , animal, day, "pyr_cycles")
    # Cycles unfiltered
    spikes=GoalFetchAnalysis.Munge.spiking.prepiso(
        spikes, l_pyr; 
        cycle=:unfilt_cycle, cells=CELLS, include_samples=false,
        matchtetrode=true)
    # Sanity check : should be more missing cycles in the filtered vs 
    #   unfiltered
    println("Missing cycles in filtered: ", sum(ismissing.(spikes.cycle)))
    println("Missing cycles in unfiltered: ", sum(ismissing.(spikes.unfilt_cycle)))
    @assert sum(ismissing.(spikes.cycle)) > sum(ismissing.(spikes.unfilt_cycle))

    spikes.isolated = Vector{Union{Bool, Missing}}(missing, nrow(spikes))

    sp = Munge.spiking.isolated!(
        subset(spikes,:animal => a->a .== animal; view=true), 
    )

end

