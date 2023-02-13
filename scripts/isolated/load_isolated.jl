begin
    using GoalFetchAnalysis
    using .Timeshift, .Plot, .Timeshift.types, .Timeshift.shiftmetrics, 
          .Field.metrics, .Plot.receptivefield, .DIutils.namedtup, 
          .Munge.isolated, .Munge.nonlocal, .Munge.spiking, .Plot.lfplot
    Filt = DI.Filt
    using .Munge.timeshift: getshift
    using .DIutils.statistic: pfunc
    filt_desc = Filt.get_filters_desc()

    using DataStructures: OrderedDict
    import DimensionalData: Between
    using ProgressMeter, DimensionalData, Infiltrator, JLD2, DataFrames,
          DataFramesMeta, StatsBase, HypothesisTests, Plots, StatsPlots

    # Parse the command line
    opt = parser()
    # Data
    data = load_iso(opt)
    lfp       = data["lfp"]; @assert(lfp isa DataFrame)
    spikes    = data["spikes"] ; @assert(lfp isa DataFrame)
    allspikes = data["allspikes"]; @assert(lfp isa DataFrame)
    tsk       = data["tsk"]; @assert(lfp isa DataFrame)
    cells     = data["cells"]; @assert(lfp isa DataFrame)
    beh       = data["beh"]; @assert(lfp isa DataFrame)
    cycles    = data["cycles"]; @assert(lfp isa DataFrame)
    
end

clab = OrderedDict(-1 => "nontask", 0 => "cue", 1=> "mem", missing=>"sleep")
Munge.nonlocal.setclab(clab)
isonames =  OrderedDict(false => :adjacent, true=>:isolated)
filt_desc = OrderedDict(:all => "> 2cm/s")
save_kws = (;pfc_rate_analy=true)
filt = Filt.get_filters()

datacut = opt["filt"]
animal  = opt["animal"]
day     = opt["day"]
tet     = opt["tet"]
Plot.setappend((;animal, day, tet))

Munge.nonlocal.setunfilteredbeh(DI.load_behavior(animal,day);
                               animal, day)

# Obtain the firing rate matrix
R = Munge.spiking.torate(allspikes, beh)

# Get PFC_units
pfc_units = @subset(cells,:area.=="PFC").unit
pfc_units = intersect(pfc_units, collect(R.dims[2]))

begin
    cycles.time = cycles.start
    @info "Registering cycles"
    DIutils.filtreg.register(cycles, spikes; on="time", transfer=["cycle"], match=:prev)
    DIutils.filtreg.register(cycles, lfp; on="time",    transfer=["cycle"], match=:prev)
    @info "Bandstopping broadraw"
    lfp = Munge.lfp.bandstop(lfp, [57, 120], [63, 750]; field=:broadraw, scale=:raw, rounds=1, order=10)
    @info "Smoothing broadraw"
    lfp=Munge.lfp.smooth(lfp, :broadraw; ker=7.5)
    @info "Plotting cycles"
    cycleplot(lfp; otherfield=[:broadraw, :filtbroadraw],
                   kws=[(;linewidth=0.25, linestyle=:dash), 
                        (;linewidth=0.5, linestyle=:dash, legend=:outerbottomright, background_color=:gray)])
    cycleplot(lfp; otherfield=[:broadraw, :smoothbroadraw, :filtbroadraw],
                   kws=[(;linewidth=0.25, linestyle=:dash), (;), (;linewidth=0.5, linestyle=:dash, legend=:outerbottomright, background_color=:gray)])
end
