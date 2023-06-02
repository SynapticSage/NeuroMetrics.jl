using GoalFetchAnalysis.Plot.lfplot
using ProgressMeter, ProfileView

animals, days = ["RY16","RY22"], [36, 21]
tetrode_sets  = [:synced, #= :best, :ca1ref =#]

global spikes, lfp, cells, opt = nothing, nothing, nothing, nothing
global cycles, ripples = nothing, nothing
global opt = nothing

iters = Iterators.product(tetrode_sets, zip(animals, days))
(tetrode_set, (animal, day)) = first(iters)

for (tetrode_set, (animal, day)) in iters
    
#Parameters
global opt = Dict(
    "animal"=>animal,
    "day"=>day,
    "load_pyr"=>false,
    "run_iso"=>false,
)
printstyled("animal = $animal, day = $day, tetrode_set = $tetrode_set",
            blink=true, color=:light_green, bold=true, reverse=true)

if isdefined(Main, :DrWatson)
    include(scriptsdir("isolated","setup_checkpoint.jl"))
else
    include("./setup_checkpoint.jl")
end
Plot.setparentfolder("isolated")

# subset(combine(groupby(cells, [:animal,:tetrode,:pyrlayer]),
#     :meanrate=>mean,
#     :meanrate=>length=>:n_cells,
#     :propbursts=>mean,
# ), :pyrlayer=>s->s.!=false)
# Time conversion factors
# lfp_convert = DI.superanimal_timeconversion("super",0)

all_tetrodes = 
if tetrode_set     == :all # all tetrodes
    tets = Dict(animal=>unique(subset(CELLS,:animal.==animal).tetrode) 
        for animal in unique(CELLS.animal))
elseif tetrode_set == :best # tetrodes with best pop phase locking across animals
    tets = Dict(
        "RY16" => [56],
        "RY22" => [29]
    )
elseif tetrode_set == :synced # tetrodes with same phase reponse across animals
    tets = Dict(
        "RY16" => [1],
        "RY22" => [29]
    )
elseif tetrode_set == :ca1ref
    tets = Dict(
        "RY16" => [:ca1ref],
        "RY22" => [:ca1ref]
    )
elseif tetrode_set == :pyr # pyramidal layer tetrodes
    DI.annotate_pyrlayer!(CELLS)
    tets = Dict(animal=>subset(CELLS, :animal=>a->a.==animal,
        :pyrlayer=>p->p.!=true, view=true).tetrode |> unique for animal in
            unique(CELLS.animal))
end

animals, days = unique(CELLS.animal), unique(CELLS.day)
ripple_props = [:ripple, :ripple_time, :ripple_phase, :ripple_phase_band,
                :ripple_amp_band]
theta_props  = [:theta, :theta_time, :theta_phase, :theta_amp]
cycle_props  = [:animal, :day, :time, :tetrode, :unit, :cycle, :unfilt_cycle]
iso_props    = [:isolated, :nearestcyc, :meancyc]
Plot.setappend("_$(animal)_$(tetrode_set)")


    println("Loading data for $animal, day $day")
       

         #  |     ,---.,---.--.--.   .,---.    |    ,---.,---.
         #  |     `---.|---   |  |   ||---'    |    |__. |---'
         #  |         ||      |  |   ||        |    |    |    
         #  `o    `---'`---'  `  `---'`        `---'`    `    
         #                                                    

        lfp = begin 
           println("Loading LFP tetrodes") 
            printstyled("tet = ", all_tetrodes[animal], "\n", blink=true)
           lfp = DI.load_lfp(animal, day, tet=all_tetrodes[animal])
           lfp.time .-= lfp_convert[animal];
           lfp
        end
        # Let's check the consistency of theta across pyr layer tetrodes
        GC.gc()
        checkranges()


         @assert !isempty(lfp) "No data!"
         # l_def = subset(lfp, :tetrode=>t->t.∈(DI.default_tetrodes[animal],),
         # view=true)
         # Get the ripple band
         fs = 1/median(diff(lfp.time[1:20_000])) 
         GoalFetchAnalysis.Munge.lfp.bandpass(groupby(lfp, :tetrode), 
                                                      150, 250; order=9, fs, 
                                                      newname="ripple")
         lfp.raw         = convert(Vector{Union{Missing, Int16}}, 
                                     round.(lfp.raw))
         lfp.ripple      = convert(Vector{Union{Missing, Float32}}, 
                                     (lfp.ripple));
         lfp.rippleamp   = convert(Vector{Union{Missing, Float32}},
                                     (lfp.rippleamp));
         lfp.ripplephase = convert(Vector{Union{Missing, Float32}},
                                     (lfp.ripplephase));

         lfp[!,:cycle] = Vector{Union{Missing, Int32}}(missing, nrow(lfp))
         @time Munge.lfp.annotate_cycles!(groupby(lfp, :tetrode), 
                                          method="peak-to-peak")
         # ISSUE: MEMORY explodes on the above line
         lfp[!,:cycle] = convert(Vector{Union{Missing, Int32}}, 
                                   lfp[!,:cycle]);
        checkranges()
        GC.gc()
        # INFO: SAVE
        lfp = DataFrame(lfp)
        lfp[!,:broadraw] = disallowmissing(lfp.broadraw);
        lfp[!,:broadraw] = convert(Vector{Float32}, lfp.broadraw);
        DI.save_lfp(transform(lfp[!,[:time, :tetrode, :cycle, :raw,
                         :phase, :broadraw, :ripple, :rippleamp, :ripplephase]],
            :time=> t-> t .+ lfp_convert[animal], renamecols=false), 
            animal, day; handlemissing=:drop, append="$tetrode_set");
        broadraw = convert(Vector{Float16}, lfp.broadraw);
        lfp = lfp[!, Not(:broadraw)];
        lfp[!,:broadraw] = broadraw;
        DIutils.pushover("...finished section I, $animal")
        checkranges()

        # ,--.     ,---.     |                                 |              
        # ,--'     `---.,---.|--- .   .,---.    ,---.,   .,---.|    ,---.,---.
        # |            ||---'|    |   ||   |    |    |   ||    |    |---'`---.
        # `--'o    `---'`---'`---'`---'|---'    `---'`---|`---'`---'`---'`---'
        #                            |             `---'                    

        # Examine
        # get second tet
        l = groupby(lfp, :tetrode)[1]
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
        cycles = Munge.lfp.get_cycle_table(groupby(lfp, :tetrode));
        cycles.cycle = convert(Vector{Union{Missing, Int32}}, cycles.cycle);
        for col in [:start, :stop, :amp_mean, :δ]
            cycles[!,col] = convert(Vector{Union{Missing, Float32}}, 
                cycles[!,col])
        end
        @assert !isempty(cycles) "No cycles!"
        @assert !all(ismissing.(cycles.amp_mean)) "No amp_mean!"
        @assert !all(ismissing.(cycles.δ)) "No δ!"
        @assert !all(ismissing.(cycles.start)) "No start!"
        @assert !all(ismissing.(cycles.stop)) "No stop!"
        # INFO: SAVE
        # Save cycles
        DI.save_cycles(begin
            transform(cycles, 
                :start=>s->s .+ lfp_convert[animal],
                :stop=>s->s  .+ lfp_convert[animal],
                renamecols=false
            )
        end, animal, day, "$tetrode_set");

        # ---------------- *
        # Visualize cycles |
        # ---------------- *
        begin
            folder, par = Plot.folder_args, Plot.parent_folder
            Plot.setparentfolder("theta");
            Plot.setfolder("cycle_cutting");
            P = []
            for g in groupby(lfp, :tetrode)
                push!(P,lfplot.cycleplot(g))
            end
            m = min(length(P),10)
            plot(P[1:m]..., layout=(m,1), size=(1200, 400*10))
            Plot.save("cycleplot_of_$(tetrode_set)_tetrodes")
            # Restore settings
            Plot.setfolder(folder...) 
            Plot.setparentfolder(par...)
        end
        checkranges()
        DIutils.pushover("Finished section II, $animal")
        current()

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
        # Ripple phase
        global spikes
        Munge.spiking.event_spikestats!(spikes, ripples; eventname="ripple",
            ampfield=:amp)
        @assert(length(unique(spikes.ripple_phase)) .> 4, 
        "Ripple phase not calculated correctly")
        @assert !all(ismissing.(spikes.ripple_phase)) "No ripple phase!"
        @assert !all(ismissing.(spikes.ripple_amp)) "No ripple amp!"
        @assert !all(ismissing.(spikes.ripple)) "No ripple!"

        # # TROUBLESHOOT:
        # # Let's figure out the number of tetrodes that match between
        # # spikes and lfp
        # sp = subset(spikes, 
        # :tetrode=>t->t .∈ (unique(l_pyr.tetrode),), view=true)
        # println("Theoretical number of pyramidal layer cells:",
        # combine(groupby(sp, :tetrode),
        # :unit => (x->length(unique(x))), renamecols=false)
        # )

        # Add ripple BAND phase (rather than phase within the ripple)
        spikes.ripple_phase_band = Vector{Union{Missing, Float32}}(missing, 
            size(spikes,1));
        spikes.ripple_amp_band   = Vector{Union{Missing, Float32}}(missing, 
            size(spikes,1));
        gl_pyr   = groupby(lfp,  :tetrode);
        g_spikes = groupby(spikes, :tetrode);
        K = intersect(keys(gl_pyr)  |>collect .|> NamedTuple,
            keys(g_spikes)|>collect           .|> NamedTuple);
        println("Number of matching keys:", length(K))

        # ISSUE: fixed badinds 5/26/23 :: rerun this
        @assert all(ripple_props .∈ (propertynames(spikes),))
        prog = Progress(length(K), 1, "Adding ripple band phase");
        k = first(K)
        for k in K
            println("Tetrode: ", k.tetrode);
            if length(K) > 1
                gs, gl = g_spikes[k], gl_pyr[k];
            else
                gs, gl = g_spikes[k], lfp;
            end
            time = gs.time
            ind = DIutils.searchsortednearest.([gl.time], time);
            try
                @assert(length(unique(ind)) > 4, "Not enough matching spikes")
            catch e
                println(e)
                @infiltrate
            end
            badinds = (gl.time[ind] - time) .> 0.03
            ind = ind[.!badinds]
            gs.ripple_phase_band[.!badinds] = gl.ripplephase[ind];
            gs.ripple_amp_band[.!badinds]   = gl.rippleamp[ind];
            next!(prog)
        end
        @assert !all(x->ismissing(x), (spikes.ripple_phase_band))
        @assert !all(x->ismissing(x), (spikes.ripple_amp_band))

        # NOTE: COMMENTED TO HERE

        # # Troubleshoot: Let's figure out what fraction of our pyramidal layer
        # #        cells
        # q=combine(groupby(sp, :tetrode),
        # :unit => (x->length(unique(x))), 
        # [:unit, :ripple_phase_band] => ((x,y)->( goodinds = .!ismissing.(y);
        #     unique(x[goodinds]) |> length)) => :n_with_ripple_phase,
        #     renamecols=false,
        # )
        # println("Theoretical number of pyramidal layer cells:", q)
        
        
        # Theta phase
        # sort!(spikes, [:unit, :time])
        # cells = subset(CELLS, :animal=>a->a.==animal, :day=>d->d.==day,
        #                 view=true)
        # DIutils.filtreg.register(cells, spikes, on="unit", transfer=["tetrode"])
        spikes.theta_amp   = Vector{Union{Missing, Float32}}(missing, size(spikes,1));
        spikes.theta_phase = Vector{Union{Missing, Float32}}(missing, size(spikes,1));
        spikes.theta       = Vector{Union{Missing, Int32}}(missing, size(spikes,1));
        Munge.spiking.event_spikestats!(groupby(spikes,:tetrode), 
                                        groupby(cycles, :tetrode); 
                                        ampfield="amp_mean",
                                        indfield="cycle",
                                        indfield_is_index=true,
                                        eventname="theta")
        @assert all(theta_props .∈ (propertynames(spikes),))
        @assert !all(x->ismissing(x), (spikes.theta_phase))
        @assert !all(x->ismissing(x), (spikes.theta_amp))
        @assert !all(x->ismissing(x), (spikes.theta))

        # q=combine(groupby(sp, :tetrode),
        # :unit => (x->length(unique(x))), 
        # [:unit, :theta_phase] => ((x,y)->( goodinds = .!ismissing.(y);
        #     unique(x[goodinds]) |> length)) => :n_with_theta_phase,
        #     renamecols=false,
        # )
        # println("Theoretical number of pyramidal layer cells:", q)
        
        # PREPARE FOR ISOLATED SPIKE ANALYSIS
        GC.gc()
        # Cycles filtered by ripples and movement
        spikes=GoalFetchAnalysis.Munge.spiking.prepiso(
            spikes, lfp; 
            cycle=:cycle,
            cells=CELLS, include_samples=false, matchtetrode=true, 
            BEH, ripples)
        # Cycles unfiltered
        spikes.unfilt_cycle = Vector{Union{Missing,Bool}}(missing, nrow(spikes))
        lfp[!,:unfilt_cycle] = lfp[!, :cycle]
        spikes=GoalFetchAnalysis.Munge.spiking.prepiso(
            spikes, lfp; 
            cycle=:unfilt_cycle, cells=CELLS, include_samples=false,
            matchtetrode=true)
        # Sanity check : should be more missing cycles in the filtered vs 
        #   unfiltered
        #   
        props = intersect(union(cycle_props, ripple_props), propertynames(spikes))
        # DI.save_spikes_taginfo(begin
        #     transform(spikes[!,props], :time=>t->t.+lfp_convert[animal], renamecols=false)
        # end, animal, day, "pyr_cycles_unfilt")
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
        # INFO: APPARENTLY runs faster if function in Main
        # @profile begin
        sp = Munge.spiking.isolated!(spikes)
        # end
        # ProfileView.view()
        # ISSUE: invesetigate long units like 134: cell has a meanrate=3.5hz)
        props = intersect(union(cycle_props, ripple_props, theta_props, iso_props),
                          propertynames(spikes))
        DI.save_spikes_taginfo(begin
            transform(spikes[!,props], 
                :time=>t->t.+lfp_convert[animal],
            renamecols=false)
        end, animal, day, "$(tetrode_set)_cycles_isolated")
        DIutils.pushover("Finished section 3, animal $animal, day $day")
end
