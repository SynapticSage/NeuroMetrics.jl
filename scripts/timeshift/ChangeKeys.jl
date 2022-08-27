using DrWatson, Revise
quickactivate(expanduser("~/Projects/goal-code"))

using PlutoUI
PlutoUI.TableOfContents(title="Caching Mains and Shuffles")

using DataFrames, DataFramesMeta
using DataStructures: OrderedDict
using KernelDensity, Distributions
using Plots, StatsPlots, Measures, Distributions
using ProgressMeter, ProgressLogging
using Combinatorics: powerset
import Base.Threads: @spawn
using ThreadSafeDicts, NaNStatistics
using Infiltrator
using TimerOutputs

using GoalFetchAnalysis 
import Timeshift
using Timeshift.dataframe: info_to_dataframe
using Field.recon_process: get_shortcutnames, inv_shortcutnames
import Load
import Filt
using Timeshift.checkpoint
using Timeshift.types: getshifts, getunits

I = load_mains()
F = load_fields()
S = load_shuffles()


# Let's change all of the main keys to one giant dataframe
Ic = copy(I)

for 🔑 ∈ keys(I)
    @info 🔑
    i = I[🔑]
    df = DataFrame()
    shifts, units = getshifts(i), getunits(i)
    for (shift,unit) in Iterators.product(shifts, units)
        if shift ∈ keys(i) .&& unit ∈ keys(i[shift])
            dfu = i[shift][unit]
            dfu[!,:shift] .= shift
            dfu[!,:unit]  .= unit[1]
            append!(df,dfu)
        end
    end
    df = hcat(df[!,[:shift, :unit]], df[!,Not([:shift,:unit])])
    I[🔑] = df
end

save_mains()

# Let's change all of the fields keys to a shiftedfields object
# (if you ever want to change the struct, you'll have to bring
# it back to a dictionary befrore this step)


# Apply shuffles to mains and fields
SIkeys = intersect(keys(I), keys(S))
