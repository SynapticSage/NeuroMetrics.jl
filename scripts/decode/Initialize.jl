using Pushover, Revise, Interact, Blink, Mux, ProgressMeter
using Statistics, NaNStatistics
using VideoIO
using GLMakie
using ColorSchemes, Colors
import ColorSchemeTools 
using DataFrames, DataFramesMeta
using Printf
using StatsPlots: @df
import StatsBase
using Infiltrator
set_theme!(theme_dark())
__revise_mode__ = :eval

using GoalFetchAnalysis
using .decode

# Debugger?
ENV["JULIA_DEBUG"] = nothing

# -----------------
# HELPER FUNCTIONS 
# -----------------
na = [CartesianIndex()]
animal, day, epoch, ca1_tetrode = "RY16", 36, 7, 5
thresh = Dict("likelihood"=>0.1, "acausal_posterior"=>0.97, "causal_posterior"=> 0.97)
typical_number_of_squares_active = quantile(1:(45*28), 1) - quantile(1:(45*28), 0.97)
transition_type, decoder_type = "empirical", "sortedspike"
variable                      = "causal_posterior"

load_from_checkpoint          = false
usevideo                      = false
remove_nonoverlap             = true # set to true if you expect multiple splits in this code
dosweep                       = false
doPrevPast                    = false
doRipplePhase                 = false
splitBehVar                   = ["egoVec_1", "egoVec_2", "egoVec_3", "egoVec_4", "egoVec_5"] # Variables to monitor with splits
vectorToWells                 = true
histVectorToWells             = true

# -----------
# Decode file locations
# -----------
dir = plotsdir("mpp_decode", "withBehVideo=$usevideo")
if !(isdir(dir))
    mkdir(dir)
end

video="/Volumes/Colliculus/RY16_experiment/actualVideos/RY16_69_0$(epoch)_CMt.1.mp4"
outputVideo = "animation.$(decoder_type)_$(transition_type)_test_0_$(variable)_$(basename(video))"
decode_file = raw.decodepath(animal, day, epoch, transition="empirical",
                             method="sortedspike", split=0,
                             type="test", speedup=20.0)
