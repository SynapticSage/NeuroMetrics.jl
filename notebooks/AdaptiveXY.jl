# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

### A Pluto.jl notebook ###
# v0.19.9
using Markdown
using InteractiveUtils
# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end


#
# !!! notebook
# 	🚀 **Adaptive receptive fields**
#
# Purpose: Test out my base adaptive field codes. Make sure the various steps (grid, occupancy, and fields) are working. Also make sure downstream shifted objects are working.
#
# ##### *TODO*
# * metrics(shifted fields)
# * shifted field plot recipes
# * radii vector support
#   * plot recipe
#   * selection via 10 increase in base sample width
# * to_dataframe
#   * fields
#   * metrics
#   * shifted fields
#
#

PlutoUI.TableOfContents(title="🚀 Adaptive RFs" )


#
# # Preamble 
# Import packages
#
#

# ╠═╡ show_logs = false
begin
	  using DrWatson
	  quickactivate(expanduser("~/Projects/goal-code"))
	  using Plots
	  using Revise
	  using DataFrames
	  using NaNStatistics
	  import ProgressLogging
	  using PlutoUI
	  using DataStructures: OrderedDict
	  using GoalFetchAnalysis
	  import Utils
	  import Timeshift
	  #import Plot
	  using Field.metrics
    using ColorSchemes
    using DataFramesMeta
	using Timeshift.types
	using Timeshift.shiftmetrics
    using StatsBase
    using DimensionalData
    using StatsPlots
    using GLM
	using ImageSegmentation, Images, LazySets, Statistics
    import Plot
	  
	  adaptive = Field.adaptive
      metrics = Field.metrics
	  WIDTHS = OrderedDict(
		  "x"=>2.5f0, "y"=>2.5f0, "currentPathLength"=>2f0,
          "currentAngle"=>Float32(2pi/80)
	  )
      filts = Filt.get_filters_precache()
	maxrad = nothing
end


prop_sel = @bind prop_str PlutoUI.Radio(["y-x","currentAngle-currentPathLength", "currentAngle","currentPathLength"], default="y-x")


#
# Loadup (spikes,beh) dataframes and transfer $(join(props,"-")) to spikes structure
#
#

beh, spikes = begin
	props = Vector{String}(split(prop_str, "-"))
	@info props
    @time spikes, beh, ripples, cells = Load.load("RY22", 21);
	@time beh, spikes  = Load.register(beh, spikes; on="time", transfer=props)
    spikes = dropmissing(spikes, props)
    beh, spikes
end;


grid_select = begin
	width_select =  @bind width  Slider(0f0:0.2f0:2f0, show_value=true, default=04f0)
	thresh_select = @bind thresh Slider(1f0:0.5f0:6f0, show_value=true, default=4f0)
	(;width_select, thresh_select)
end


#
#     widths = $widths
#     
#

radiusinc, ylim, aspect_ratio = if prop_sel == "y-x"
	0.1f0, nothing, 1
elseif prop_sel == "currentAngle-currentPathLength"
    0.05f0, (0, 100), 1/18
else
	0.05f0, nothing, 1
end

widths = [WIDTHS[prop]*width for prop in props]


#
# # Δt = 0 only
#
# ## 🌐 Grid
# Let's try making a grid object
#
#

widths


# ╠═╡ show_logs = false
begin
    boundary = prop_str == "y-x" ? nothing : Dict("currentPathLength"=>(0,100))
    @time G = adaptive.get_grid(beh, props; widths, thresh, maxrad=prop_str == "y-x" ? 10f0 : nothing, radiusinc, boundary);
end


#
# Implied linear width of `maxrad[$(nanmaximum(G.radii))]` => $(nanmaximum(G.radii) * sqrt(2)) 
#
#

grid_select


#plot(G; title="radii\nresolution=$(size(G.grid))")


#if ndims(G.radii) == 2
#	heatmap([collect(x) for x in G.centers]..., (G.radii .=== NaN32)'; title="nan locations")
#end


G.centers


G.radii


unique(G.radii)


#
# ## 💠 Occupancy
# And an occupancy object
#
#

O = @time adaptive.get_occupancy(beh, G);


#plot(O, clim=(0,0.01), ylims=ylim)




# ╠═╡ disabled = true
#=╠═╡
O
  ╠═╡ =#


#
# ## Multiunit adaptive field
# grab all spikes as single multiunit field
#
#

# Test field abilities
@time multiunit = @time adaptive.get_adaptivefield(spikes, G, O);
# @benchmark adaptive.get_adaptivefield(spikes, G, O);


grid_select


#plot(multiunit)


#
# ## 🔥 Field dict of all cells
#
#

begin
	@time units = adaptive.yartsev(spikes, G, O; widths=width, thresh, 
	                               filters=filts[:all]);
    unitsγG = units
end;
grid_select


#
# ### Visualize fields @ Δt=0
U = units[(;unit=unit)];


revise(Plot.receptivefield)


begin
	bps = OrderedDict(k=>v[:bitsperspike] for (k,v) in units)
	sort!(bps, by=k->bps[k], rev=true)
	units_ordered = [x[1] for x in keys(bps)]
end


# ╠═╡ disabled = true
#=╠═╡
Memoization.empty_all_caches!()
  ╠═╡ =#


# ╠═╡ disabled = true
#=╠═╡
Plot.setfolder("goal", "pathlength")
  ╠═╡ =#


unit_select = @bind unit PlutoUI.Slider(units_ordered, default=31, show_value=true)


plot(U)


# ╠═╡ disabled = true
#=╠═╡
Plot.save((;desc="pathlength",unit))
  ╠═╡ =#


ylims === nothing


μ_firing = begin
    Q = units[(;unit=unit)]
    nansum(reshape(Q.occ.prob, size(Q.occ.count)) .* Q.rate)
end


field = units[(;unit=unit)]

