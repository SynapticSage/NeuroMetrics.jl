### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ fc614ab8-00cb-11ed-0f62-f751ef056b39
begin
	
	using DrWatson, Revise
	using PlutoUI
	quickactivate(expanduser("~/Projects/goal-code"))
	using DataFrames, DataFramesMeta
	using DataStructures: OrderedDict
	using KernelDensity, Distributions
	using Plots, StatsPlots, Measures, Distributions
	using ProgressMeter, ProgressLogging
	using Combinatorics: powerset
	import Base.Threads: @spawn
	using ThreadSafeDicts, NaNStatistics
	

	using GoalFetchAnalysis
	import Timeshift
	using Timeshift.dataframe: info_to_dataframe
	using Field.recon_process: get_shortcutnames, inv_shortcutnames
	import Load
	import Filt
	filts = Filt.get_filters()

	function ingredients(path::String)
		# this is from the Julia source code (evalfile in base/loading.jl)
		# but with the modification that it returns the module instead of the last object
		name = Symbol(basename(path))
		m = Module(name)
		Core.eval(m,
	        Expr(:toplevel,
	             :(eval(x) = $(Expr(:core, :eval))($name, x)),
	             :(include(x) = $(Expr(:top, :include))($name, x)),
	             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
	             :(include($path))))
		m
	end
    sets = ingredients(scriptsdir("timeshift", "TimeShift_setsOfInterest.jl"))
	prop_set = sets.marginals_superhighprior_shuffle


    function get_key(;shifts, kws...)
        (;kws..., grid=:adaptive,
		first=first(shifts), last=last(shifts), 
				   step=Float64(shifts.step)) 
    end
    function keymessage(I::AbstractDict, key)
        docontinue=false
        if key ∈ keys(I)
            if I[key] isa Task && !(istaskfailed(I[key]))
                #@info "task key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                docontinue=true
            elseif I[key] isa Task && istaskfailed(I[key])
                "key=$key already exists, but failed...redo!"
            else
                #@info "key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                docontinue=true
            end
        end
        if key ∉ keys(I)
            #@info "key=$key ∉ keys, ...creating..."
        end
		docontinue
     end
     PROPS = ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
     IDEALSIZE = Dict(key => (key=="stopWell" ? 5 : 40) for key in PROPS)
	 shifts=2:0.05:2
end

# ╔═╡ c99b4903-464d-44fb-b5e5-d7724b25afea
PlutoUI.TableOfContents(title="Caching Mains and Shuffles")

# ╔═╡ 823b1bff-d922-4c2b-8a50-179af24094bd
begin
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
	_, spikes = Load.register(beh, spikes; transfer=["velVec"], on="time")
end;

# ╔═╡ 7ff3160c-fbf4-4740-8b24-5e3a2f130f70
begin
	datacuts = collect(keys(filts))
	(;prop_set, datacuts)
end

# ╔═╡ a0574649-00a6-4883-9a87-0dac3fc5f8a5
md"""
# Main
"""

# ╔═╡ caaff430-fa09-4bbd-b106-95e36d9203a3
md"## Load checkpoint"

# ╔═╡ 90dfe32b-0930-4067-8936-6f1e1e922a35
begin
    if isfile(Timeshift.mainspath())
        I = Timeshift.load_mains()
    else
        I = OrderedDict()
    end
	keys(I)
end

# ╔═╡ c820cf54-fa0a-4112-8e76-8f76839b7a49
md"## Cache results"

# ╔═╡ 673b09d2-5dd4-4b6c-897e-2fc43f04ab8f
# ╠═╡ show_logs = false
begin
    @progress "Datacut iteration" for datacut ∈ datacuts
        finished_batch = false
        @progress "Props" for props ∈ prop_set
            marginal = get_shortcutnames(props)
            key = get_key(;marginal, datacut, shifts)
			#@info key
    		if keymessage(I, key); continue; end
            I[key] = Timeshift.shifted_fields(beh, spikes, shifts, props; widths=2.50f0)
            finished_batch = true
        end
        if finished_batch
            Timeshift.save_mains(I)
        end
    end

end

# ╔═╡ 08b1ac9b-58e8-41de-994f-a05609df3b2c
keys(I)

# ╔═╡ cfb3aaab-2166-43ce-9fdf-586b98fe8272
md"""
## Test a single key
Works for a single key?"""

# ╔═╡ 135856f2-6c6b-4bc4-9e7a-ca678d5e729d
# ╠═╡ disabled = true
#=╠═╡
datacut, props = first(datacuts), first(prop_set)
  ╠═╡ =#

# ╔═╡ 2f11999f-5d6e-4807-8bce-49791e7a0211
# ╠═╡ disabled = true
#=╠═╡
marginal=get_shortcutnames(props)
  ╠═╡ =#

# ╔═╡ ae54aa7a-3afb-43c9-b083-0908b1f02d18
#=╠═╡
key = get_key(;marginal, datacut, shifts)
  ╠═╡ =#

# ╔═╡ dde2b892-eae4-4dfd-8735-b533e8a5ae68
#=╠═╡
Timeshift.shifted_fields(beh, spikes, shifts, props; widths=2.50f0)
  ╠═╡ =#

# ╔═╡ 835d1acf-95f6-47c1-9bdd-0b0a75034353
md"""
# Shuffles

## Load checkpoint
"""

# ╔═╡ 29497f7b-795e-433f-b772-72191f52dc24
begin
    if isfile(Timeshift.mainspath())
        S = Timeshift.load_shuffles()
    else
        S = OrderedDict()
    end
	keys(S)
end

# ╔═╡ 2696b1df-e52c-497c-b62f-a0932da6c8a4


# ╔═╡ Cell order:
# ╟─c99b4903-464d-44fb-b5e5-d7724b25afea
# ╟─fc614ab8-00cb-11ed-0f62-f751ef056b39
# ╟─823b1bff-d922-4c2b-8a50-179af24094bd
# ╟─7ff3160c-fbf4-4740-8b24-5e3a2f130f70
# ╟─a0574649-00a6-4883-9a87-0dac3fc5f8a5
# ╟─caaff430-fa09-4bbd-b106-95e36d9203a3
# ╠═90dfe32b-0930-4067-8936-6f1e1e922a35
# ╟─c820cf54-fa0a-4112-8e76-8f76839b7a49
# ╠═673b09d2-5dd4-4b6c-897e-2fc43f04ab8f
# ╠═08b1ac9b-58e8-41de-994f-a05609df3b2c
# ╟─cfb3aaab-2166-43ce-9fdf-586b98fe8272
# ╠═135856f2-6c6b-4bc4-9e7a-ca678d5e729d
# ╠═2f11999f-5d6e-4807-8bce-49791e7a0211
# ╠═ae54aa7a-3afb-43c9-b083-0908b1f02d18
# ╠═dde2b892-eae4-4dfd-8735-b533e8a5ae68
# ╟─835d1acf-95f6-47c1-9bdd-0b0a75034353
# ╠═29497f7b-795e-433f-b772-72191f52dc24
# ╠═2696b1df-e52c-497c-b62f-a0932da6c8a4
