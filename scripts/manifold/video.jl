using DrWatson
quickactivate(expanduser("~/Projects/goal-code"))
using GoalFetchAnalysis, Munge.causal, Munge.manifold
import Load, Munge

using Revise, ProgressMeter, CategoricalArrays, Statistics, NaNStatistics, VideoIO, GLMakie, ColorSchemes, Colors, DataFrames, DataFramesMeta, Printf, Infiltrator, ImageFiltering
using StatsPlots: @df
import StatsBase, ColorSchemeTools 


## ----------
## PARAMETERS
## ----------
animal, day, filt, N = "RY16", 36, :all, 100
areas = (:ca1,:pfc)
distance = :many
feature_engineer = :many # many | nothing
esttype = :binned
est, params = get_est_preset(esttype, horizon=1:60, thread=true, binning=7, window=1.25)

# Load all requisite vars
manifold.load_manis_workspace(Main, animal, day; filt, 
      areas, distance, feature_engineer, 
      N)
spikes, beh, ripples, cells  = Load.load(animal, day)
storage = load_alltimes_savefile(animal, day, N; params)
video             = Load.video.load_videocollection(animal, day)
TS = [Load.video.load_videots(animal, day, i) for i in 1:8]
exampframe = video[0.0]
xax, yax = collect.(exampframe.dims)
if (:index ∉ propertynames(beh))
    beh[!,:index] = 1:size(beh,1)
end

# -------------------------
# Filter by some conditions
# -------------------------
# Behavior
beh = @subset(beh, :epoch .== 2)
T = size(beh, 1)

# Filter manifold
ems = @subset(em, 
              :metric .== "CityBlock",
              :feature .== :raw,
              :dim .== 3)
