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
using DimensionalData, Distributions

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
    if !(typeof(i) <: AbstractDataFrame)
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
    end
    I[🔑] = df
end

save_mains(I)

# Let's change all of the main keys to one giant dataframe
Ic = copy(S)

for 🔑 ∈ keys(S)
    @info 🔑
    s = S[🔑]
    s = OrderedDict{keytype(s)}{Union{DataFrame,valtype(s)}}(s)
    for shuf in keys(s)
        i = s[shuf]
        if !(typeof(i) <: AbstractDataFrame)
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
        end
        s[shuf] = df
    end
    S[🔑] = s
end

save_shuffles(S; overwrite=true)

# Let's change all of the fields keys to a shiftedfields object
# (if you ever want to change the struct, you'll have to bring
# it back to a dictionary befrore this step)

Ic = copy(F)

for 🔑 ∈ keys(F)
    @info 🔑
    f = F[🔑]
    shifts, units = getshifts(f), getunits(f)
    if !(typeof(f) <: ShiftedFields)
        f = ShiftedFields(f)
    end
    F[🔑] = f
end

save_fields(F; overwrite=true)

# Apply shuffles to mains and fields
SIkeys = intersect(keys(I), keys(S))

🔑 = collect(SIkeys)[2]
i, s    = I[🔑], S[🔑]


@time s = Table.to_dataframe(s)

s = unstack(s[!, hasproperty(s, :dim_1) ? Not(:dim_1) : Colon()], :metric, :value)
i = unstack(i[!,Not(:dim_1)], :metric, :value)
Table.vec_arrayofarrays!(s)
Table.vec_arrayofarrays!(i)

# Create gaussian to sample from
shifts, units = unique(i.shift), unique(i.unit)
dists = DimArray(
         Array{Distribution,2}(undef, length(shifts), length(units)), 
         (Dim{:shift}(shifts), Dim{:unit}(units))
        );
for (shift,unit) in Iterators.product(shifts,units)
    
    ss = @subset(s, :shift .== shift, :unit .== unit)
    ii = @subset(i, :shift .== shift, :unit .== unit)
    vals = ss.bitsperspike
    μ, σ = mean(vals), std(vals)
    @infiltrate μ > 0
    D = Distributions.Gaussian(μ, σ)
    dists[unit=At(unit), shift=At(shift)] = D

end

