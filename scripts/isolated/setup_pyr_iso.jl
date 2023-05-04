@time using GoalFetchAnalysis
using DataFrames, Statistics
using Plots
animal, day = "super_clean", 0
spikes, beh, cells, ripples = DI.load(animal, day, 
    data_source=["spikes","behavior", "cells", "ripples"])
spikes.cycle = Vector{Union{Float32, Missing}}(missing, size(spikes,1)) 
DI.annotate_pyrlayer!(cells)
Plot.setparentfolder("isolated")

subset(combine(groupby(cells, [:animal,:tetrode,:pyrlayer]),
    :meanrate=>mean,
    :meanrate=>length=>:n_cells,
    :propbursts=>mean,
), :pyrlayer=>s->s.!=false)


# Time conversion factors
time_factors = DI.superanimal_timeconversion("super",0)
animals, days = unique(cells.animal), unique(cells.day)
animal, day = zip(animals, days)|>first
DIutils.pushover("Finished loading non-lfp")
animal = animals[1]

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
    # Let's check the consistency of theta across pyr layer tetrodes
    all_tetrodes = pyr_tetrodes = subset(cells, :pyrlayer=>s->s.!=false,
        :animal=>a->a.==animal).tetrode |> unique
    all_tetrodes = [all_tetrodes..., DI.default_tetrodes[animal],
                    DI.ca1ref_tetrodes[animal]]
    @time lfp  = subset(lfp, :tetrode=>t->t.∈(all_tetrodes,), view=true);
    ulfp = unstack(lfp, :time, :tetrode, :raw);
    sort(Int.(unique(lfp.tetrode)))
    if !isempty(ulfp)
        DIutils.pushover("Finished loading lfp for $animal")
    else
        DIutils.pushover("No lfp for $animal")
    end

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
    current()

    l_pyr = subset(lfp, :tetrode=>t->t.∈(pyr_tetrodes,), view=true);
    # l_def = subset(lfp, :tetrode=>t->t.∈(DI.default_tetrodes[animal],),
    # view=true)
    
    # setup cycles
    @time Munge.lfp.annotate_cycles!(groupby(l_pyr, :tetrode), 
                                     method="peak-to-peak")
    l_pyr[!,:cycle] = convert(Vector{Union{Missing, Int32}}, l_pyr[!,cycle]);
    GC.gc()
    DI.save_lfp(transform(l_pyr[!,[:time, :tetrode, :cycle]],
        :time=> t-> t .+ time_factors[animal]), animal, day; append="pyr")
    DIutils.pushover("Finished annotating cycles for $animal")

    cycles = Munge.lfp.get_cycle_table(groupby(l_pyr, :tetrode))
    cycles.cycle = convert(Vector{Union{Missing, Int32}}, cycles.cycle)
    for col in [:start, :stop, :amp_mean, :δ]
        cycles[!,col] = convert(Vector{Union{Missing, Float32}}, cycles[!,col])
    end


    DI.save_cycles(cycles, animal, day, "pyr")

    # α Setup spikes
    Munge.spiking.prepiso!(
        subset(spikes,:animal => a->a .== animal; view=true), 
        l_pyr; cells=cells, include_samples=false, matchtetrode=true, beh,
        ripples)

end

# Below I have a short code blurb that can reload the elements in the α labeled 
# comment above
begin
    spikes, beh, cells, ripples = DI.load(animal, day, 
        data_source=["spikes","behavior", "cells", "ripples"])
    l_pyr = DI.load_lfp(animal, day; append="pyr")
    DI.cyclepath(animal, day, "pyr")
    cycles = DI.load_cycles(animal, day, "pyr")
end
