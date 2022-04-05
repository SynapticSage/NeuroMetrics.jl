# LFP scratchpad
# for playing around with techinques for processing theta cycles

@time lfpd, lfp = begin
    lfp = combine(groupby(lfp, :tetrode), raw.lfp.annotate_cycles)
    lfpd = combine(groupby(lfp, :tetrode),
                   df->raw.downsample(df;dfactor=15))#little more memory hungry than for-loop version
    #lfpd = combine(groupby(lfpd, :tetrode), raw.lfp.annotate_cycles)
    tetrode, lfpd = raw.register(tetrode, lfpd; transfer=["area"], on="tetrode");
    #lfp_avg = combine(groupby(lfpd, :area), raw.lfp.mean_lfp)
    #@time lfp_wavg = combine(groupby(lfpd, :area), raw.lfp.weighted_lfp)
    #lfp_avg = combine(groupby(lfp_avg, :area), raw.lfp.annotate_cycles)
    lfpd, lfp
end
cycle_max = combine(groupby(lfpd, [:tetrode, :area]), :cycle=>maximum)
ca1, pfc = groupby(lfpd, :area)
uca1 = raw.lfp.unstack_tetrode(ca1)
upfc = dropmissing(raw.lfp.unstack_tetrode(pfc))

# Separate out groups of correlated theta phase C>1 and average those (until this is done, use tet=1 to cut theta cycles)
import NaNStatistics, Missings
CA1_c  = NaNStatistics.nancor(Matrix(uca1[:, 2:end]))
good_ca1_tets = parse.(Int, names(uca1)[2:end][CA1_c[1,:] .> 0])
filt = in.(ca1.tetrode, [good_ca1_tets])
ca1_avg = raw.lfp.mean_lfp(ca1[filt,:])

# Butterworth and hilbert the averaged raw (averaging introduces higher
# freequency changes)
filt = DSP.analogfilter(DSP.Bandpass(6, 12, fs=1/median(diff(ca1_avg.time))),
                 DSP.Butterworth(5))
ca1_avg.raw, ca1_avg.phase  =  Float64.(ca1_avg.raw), Float64.(ca1_avg.phase)
#x = DSP.filtfilt(filt.p, ca1_avg.phase)
hilb = DSP.hilbert(ca1_avg.raw)
ca1_avg.amp, ca1_avg.phase = abs.(hilb), angle.(hilb)

# Find higher variance tetrodes
ca1_avg = raw.lfp.gauss_lfp(ca1_avg)
filt = in.(ca1.tetrode, [good_ca1_tets])
variances = combine(groupby(ca1[filt,:],:tetrode),:raw=>var)
ranked_variances = sort(variances,:raw_var,rev=true)


# Make cycle table for high var tetrodes
function get_cycle_table(lfp)
    @assert "cycle" in names(lfp)
    tab = table.get_periods(lfp, "cycle")
    # TODO will eventually add more stats before this function exits
    return tab
end
cycles = get_cycle_table(filter(:tetrode=> t->t==5, lfp))

## PLOTTING ###
# Checking the lfp
using StatsPlots
@df cycle_max Plots.scatter(:area, :cycle_maximum, group=:area)
@df cycle_max Plots.histogram(:cycle_maximum)
@df lfpd Plots.plot(:time, :cycle, group=:tetrode)
@df uca1[1:10:10000,:] corrplot(cols(2:10))


Plots.heatmap(CA1_c, c=:vik, clims=(-1,1), xticks=(1:size(C,2),
                                                   names(uca1)[2:end]),
              xrotation=90, yticks=(1:size(C,2), names(uca1)[2:end]),
              yrotation=0, title="CA1 Theta Phase Correlation")
Plots.savefig(plotsdir("theta","heatmap_correlation_thetaphase_ca1.svg"))

C = Statistics.cor(Matrix(upfc[:, 2:end]))
Plots.heatmap(C, c=:vik, clims=(-1,1), xticks=(1:size(C,2),
                                               names(uca1)[2:end]),
xrotation=90, yticks=(1:size(C,2), names(uca1)[2:end]), yrotation=0,
title="Theta Phase Correlation")
Plots.savefig(plotsdir("theta","heatmap_correlation_thetaphase_pfc.svg"))
