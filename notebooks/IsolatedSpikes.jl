### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 17d367e4-2823-11ed-10c2-ed3e82826650
begin
    using DrWatson
    quickactivate(expanduser("~/Projects/goal-code"))
end

# ╔═╡ 65d9113e-254a-4f2a-bea2-b26bb48cdd1f
begin
	 using GoalFetchAnaylsis
		using Timeshift
    using Timeshift.types
    using Timeshift.shiftmetrics
    using Field.metrics
    using Plot
    using Plot.receptivefield
    using Utils.namedtup
    using DataStructures: OrderedDict
    using DimensionalData
    import DimensionalData: Between
    using ProgressMeter
    import Plot
    using Munge.spiking
    using Filt
    using DataFrames, DataFramesMeta
    using Statistics, NaNStatistics, StatsBase, StatsPlots, HypothesisTests, GLM
    using Plots
    #using ElectronDisplay
    using LazySets
    using Munge.timeshift: getshift
    using Utils.stats: pfunc
    clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
    isonames =  OrderedDict(false => :adjacent, true=>:isolated)
    filt_desc = OrderedDict(:all => "> 2cm/s")
    save_kws = (;pfc_rate_analy=true)
    filt = Filt.get_filters()
end

# ╔═╡ 17da6d36-2823-11ed-2b22-9532d8f2b9bd

begin
    # Acquire data
    @time spikes, beh, cells = Load.load("RY16", 36, data_source=["spikes","behavior", "cells"])
    beh, spikes = Utils.filtreg.filterAndRegister(beh, spikes, on="time", transfer=["x","y","cuemem"], filters=filt[:all], filter_skipmissingcols=true)
    allspikes = copy(spikes)
    beh2 = Load.load_behavior("RY16",36)

    # Acquire LFP and isolated spikes
    lfp = Load.load_lfp("RY16", 36, tet=5);
    lfp.time = lfp.time .- Load.min_time_records[1]
    lfp = Munge.lfp.annotate_cycles(lfp)
    #sp = @subset(spikes, :tetrode .== 6);
    Munge.spiking.isolated(spikes, lfp)
end


# ╔═╡ b868cf8d-96db-49ab-b678-caeda3878583
begin
	Plot.setfolder("nonlocality","MUA and isolated MUA")
	kws=(;legend_position=Symbol("outerbottomright"))
	
end

# ╔═╡ 4379b31d-cb3b-43fb-b983-30c60c223e47
begin
	
	@df @subset(iso_sum,:area.=="CA1") bar(:cmlab, :timespent, ylabel="Time spent", group=:cuemem; kws...)
	Plot.save((;desc="time spent"));
end

# ╔═╡ 43e913fd-6ab3-478a-96a8-de3a2bf2bf2a

@df iso_sum bar(:cuearea, :events_per_time, ylabel="MUA events per second\n$(filt_desc[:all])", group=:cuemem; kws...)
Plot.save((;desc="MUA per second"));


# ╔═╡ 5a9c4f7b-29d1-4083-ada9-7b0a3987c038

@df iso_sum bar(:cuearea, :isolated_mean, ylabel="Isolated Spikes (sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:cmlab; kws...)
Plot.save((;desc="fraction of isolated spikes"));

# ╔═╡ 6483cf7b-e633-4e68-ba84-f55096777407

@df iso_sum bar(:cuearea, :isolated_events_per_time, ylabel="Isolated MUA × sec⁻1\n(sign of CA1-PFC interaction)\n$(filt_desc[:all])", group=:cuemem; kws...)
Plot.save((;desc="isolated spikes per second"));

# ╔═╡ d975603c-8461-4ade-8357-bf154640d439


# ╔═╡ f62a865f-87e7-4ee1-85fc-e18fc4bdaf59


# ╔═╡ 2a837293-3dc2-4306-9c5f-2c16d75c83d6


# ╔═╡ 8ac53f41-ef29-488f-9338-58ed05e879f4


# ╔═╡ 6c4db7c0-24ca-44f8-969c-ace3db3466e5


# ╔═╡ 22ec7458-eb52-4563-846a-eac204709f9c


# ╔═╡ a18d17ba-e4a9-4810-a2fe-5e54cc78850a


# ╔═╡ 00ac21e3-f1b3-45d0-9319-57eecf51a252


# ╔═╡ 352449f4-47f1-4ee3-9b80-c0451723662d


# ╔═╡ 4931756b-dd5f-4b1b-92bd-60bc9b4a9590


# ╔═╡ 72c01dab-8259-4ca4-91ad-69198be9e553


# ╔═╡ bdc8e591-0fa1-463f-88ae-07835b040893


# ╔═╡ 65142c62-d204-4140-990f-0269b7e1880a


# ╔═╡ cff64793-7104-4c5e-b51a-48faec4712d1


# ╔═╡ fe56e337-1b1f-459e-b394-4f7d4ae5738d


# ╔═╡ 2e1481b0-6ef7-4703-9eca-1f0d5a3b4252


# ╔═╡ 784bb992-db1f-4307-a190-3b5136d96788
begin
# Acqruire events per time as events  / time spent
iso_sum.events_per_time = iso_sum.x1 ./ (iso_sum.timespent)
iso_sum.cuearea = iso_sum.area .* "\n" .* getindex.([clab], iso_sum.cuemem)
iso_sum.cmlab = getindex.([clab], iso_sum.cuemem)
iso_sum.isolated_events_per_time = iso_sum.isolated_mean .* iso_sum.events_per_time

ord = Dict("nontask"=>1,"cue"=>2,"mem"=>3)
iso_sum = sort(iso_sum, [DataFrames.order(:cmlab, by=x->ord[x]),:cuearea])

pfc_units = @subset(cells,:area.=="PFC").unit
R = Munge.spiking.torate(allspikes, beh)
pfc_units = intersect(pfc_units, collect(R.dims[2]))

isolated = last(groupby(subset(spikes, :isolated=>x->(!).(isnan.(x))) ,
                                        :isolated))
@assert all(isolated.isolated .== true)

end

# ╔═╡ 17da6d70-2823-11ed-30aa-2144503b97cd
begin
# Split by isolated spikes and discover the fraction of isolated spikes
iso_sum = combine(groupby(spikes, [:area, :cuemem]), :isolated => mean, (x->nrow(x)))

# Calculate time animal spends in each cuemem segment
task_pers = Table.get_periods(beh2, [:traj, :cuemem], timefract= :velVec => x->abs(x) > 2)

# Total that time and register that column to the isolation summary
task_pers = combine(groupby(task_pers, [:cuemem]), [:δ,:frac] => ((x,y)->sum(x.*y)) => :timespent)
Utils.filtreg.register(task_pers, iso_sum, on="cuemem", transfer=["timespent"])

end


# ╔═╡ Cell order:
# ╠═17d367e4-2823-11ed-10c2-ed3e82826650
# ╠═65d9113e-254a-4f2a-bea2-b26bb48cdd1f
# ╠═17da6d36-2823-11ed-2b22-9532d8f2b9bd
# ╠═17da6d70-2823-11ed-30aa-2144503b97cd
# ╠═784bb992-db1f-4307-a190-3b5136d96788
# ╠═b868cf8d-96db-49ab-b678-caeda3878583
# ╠═4379b31d-cb3b-43fb-b983-30c60c223e47
# ╠═43e913fd-6ab3-478a-96a8-de3a2bf2bf2a
# ╠═5a9c4f7b-29d1-4083-ada9-7b0a3987c038
# ╠═6483cf7b-e633-4e68-ba84-f55096777407
# ╠═d975603c-8461-4ade-8357-bf154640d439
# ╠═f62a865f-87e7-4ee1-85fc-e18fc4bdaf59
# ╠═2a837293-3dc2-4306-9c5f-2c16d75c83d6
# ╠═8ac53f41-ef29-488f-9338-58ed05e879f4
# ╠═6c4db7c0-24ca-44f8-969c-ace3db3466e5
# ╠═22ec7458-eb52-4563-846a-eac204709f9c
# ╠═a18d17ba-e4a9-4810-a2fe-5e54cc78850a
# ╠═00ac21e3-f1b3-45d0-9319-57eecf51a252
# ╠═352449f4-47f1-4ee3-9b80-c0451723662d
# ╠═4931756b-dd5f-4b1b-92bd-60bc9b4a9590
# ╠═72c01dab-8259-4ca4-91ad-69198be9e553
# ╠═bdc8e591-0fa1-463f-88ae-07835b040893
# ╠═65142c62-d204-4140-990f-0269b7e1880a
# ╠═cff64793-7104-4c5e-b51a-48faec4712d1
# ╠═fe56e337-1b1f-459e-b394-4f7d4ae5738d
# ╠═2e1481b0-6ef7-4703-9eca-1f0d5a3b4252
