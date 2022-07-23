### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ c99b4903-464d-44fb-b5e5-d7724b25afea
begin
	using DrWatson, Revise
	quickactivate(expanduser("~/Projects/goal-code"))
	using PlutoUI
	PlutoUI.TableOfContents(title="Caching Mains and Shuffles")
end

# ╔═╡ fc614ab8-00cb-11ed-0f62-f751ef056b39
# ╠═╡ show_logs = false
begin

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
	filts = Filt.get_filters_precache()

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
        @info key
        docontinue=false
        if Utils.namedtup.orderlessmatch(key, keys(I))
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
	 shifts=-2:0.05:2
	 widths = 5.0f0
	 thresh = 1.5f0

	shuffle_type = :dotson
	datacuts = collect(keys(filts))

end;

# ╔═╡ 0f48f044-fd14-4b15-bf0c-b39c1843f9db
md"""
### Caching Mains and Shuffles

Purpose: This notebook functions to cache shifted fields and shifted field shuffles under certain settings.
"""

# ╔═╡ 823b1bff-d922-4c2b-8a50-179af24094bd
# ╠═╡ show_logs = false
begin
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
	_, spikes = Load.register(beh, spikes; transfer=["velVec"], on="time")
	if shuffle_type == :dotson
	    nbins = 50
	    Munge.behavior.annotate_relative_xtime!(beh)
	    beh.trajreltime_bin = floor.(beh.trajreltime * (nbins-1))
	    _, spikes = Load.register(beh, spikes;
	                             transfer=["trajreltime","trajreltime_bin"],
	                             on="time")
	end
end;

# ╔═╡ f6add475-5405-42d1-a71c-e309aaf0121e
md"""
### Parameter description

Running shifting procedure w/:
"""

# ╔═╡ 0bf4d984-bd18-4dd2-905b-04cb5d682556
begin
	tab="""
\\begin{aligned}
&\\begin{array}{cc}
\\hline \\hline \\text { Param } & \\text { Values } \\\\
\\hline 
\\text {shifts} & $shifts \\\\
\\text {widths} & $widths \\\\
\\text {shuffle type} & \\text {$(String(shuffle_type))} \\\\
\\text {prop set} & \\text {$(Symbol(prop_set))}\\\\
\\text {seconds of sample req} & \\text {$(Symbol(thresh))}\\\\
\\hline
\\end{array}
\\end{aligned}
""";
	md"$(Markdown.LaTeX(tab))"
end

# ╔═╡ 698f5a55-2daf-4ca5-9f84-a889a491acc5
md"""
### Filtration conditions for data cuts
"""

# ╔═╡ 690df6ec-71ec-46cd-8699-960217bd4f06
filts

# ╔═╡ a0574649-00a6-4883-9a87-0dac3fc5f8a5
md"""
# Main
"""

# ╔═╡ caaff430-fa09-4bbd-b106-95e36d9203a3
md"""## Load checkpoint
Obtain mains checkpoint data
"""

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
begin
    @progress "Datacut iteration" for datacut ∈ datacuts
        finished_batch = false
        for props ∈ prop_set
            marginal = get_shortcutnames(props)
            key = get_key(;marginal, datacut, shifts, widths, thresh)
			#@info key
			filt = filts[datacut]
			if filt === nothing
				continue
			end
    		#if keymessage(I, key); continue; end
            I[key] = Timeshift.shifted_fields(beh, spikes, shifts, props; widths, filters=filt, thresh)
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
datacut, props = first(datacuts), first(prop_set)

# ╔═╡ 2f11999f-5d6e-4807-8bce-49791e7a0211
marginal=get_shortcutnames(props)

# ╔═╡ ae54aa7a-3afb-43c9-b083-0908b1f02d18
key = get_key(;marginal, datacut, shifts, widths, thresh)

# ╔═╡ d92dc54e-2ba6-4fa1-bd43-7d06c5d9cb6f
single_filt = filts[datacut]

# ╔═╡ 93a9beb2-084d-4fe7-937c-39d74740cade
md"""
### Shuffle
"""

# ╔═╡ 5941ed8d-62c0-4d11-b897-30fc24f25b78
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
Timeshift.shuffle.shifted_field_shuffles(beh, spikes, shifts, props; 
fieldpreset=:yartsev, shufflepreset=shuffle_type, nShuffle=3,
widths=2.50f0)
  ╠═╡ =#

# ╔═╡ f8f1c3e8-ed99-4a37-8faa-82ec1f946898
md"""
### Main
"""

# ╔═╡ dde2b892-eae4-4dfd-8735-b533e8a5ae68
# ╠═╡ disabled = true
#=╠═╡
@time Timeshift.shifted_fields(beh, spikes, shifts, props; widths=2.50f0, filters=single_filt, thresh);
  ╠═╡ =#

# ╔═╡ e4acf92a-3242-4579-a046-95649f835c36
;

# ╔═╡ 835d1acf-95f6-47c1-9bdd-0b0a75034353
md"""
# Shuffles

## Load checkpoint
"""

# ╔═╡ 29497f7b-795e-433f-b772-72191f52dc24
# ╠═╡ disabled = true
#=╠═╡
begin
    if isfile(Timeshift.shufflespath())
        S = Timeshift.load_shuffles()
    else
        S = OrderedDict()
    end
	keys(S)
end
  ╠═╡ =#

# ╔═╡ 2696b1df-e52c-497c-b62f-a0932da6c8a4
#=╠═╡
begin
    @showprogress "Datacut iteration" for datacut ∈ datacuts
        finished_batch = false
        @showprogress "Props" for props ∈ prop_set
            marginal = get_shortcutnames(props)
            key = get_key(;marginal, datacut, shifts, widths, thresh)
            filt = filts[datacut]
    		if keymessage(S, key); continue; end
            @time S[key] = Timeshift.shuffle.shifted_field_shuffles(beh, spikes, shifts, props; fieldpreset=:yartsev, shufflepreset=shuffle_type, nShuffle=100, widths, thresh, exfiltrateAfter=10, shiftbeh=false, filters=filt)
            finished_batch = true
        end
        if finished_batch
            Timeshift.save_shuffles(S)
        end
    end
end
  ╠═╡ =#

# ╔═╡ eb7eb4ca-6088-487c-9017-6b7988188c20


# ╔═╡ Cell order:
# ╟─0f48f044-fd14-4b15-bf0c-b39c1843f9db
# ╟─c99b4903-464d-44fb-b5e5-d7724b25afea
# ╟─fc614ab8-00cb-11ed-0f62-f751ef056b39
# ╟─823b1bff-d922-4c2b-8a50-179af24094bd
# ╟─f6add475-5405-42d1-a71c-e309aaf0121e
# ╟─0bf4d984-bd18-4dd2-905b-04cb5d682556
# ╟─698f5a55-2daf-4ca5-9f84-a889a491acc5
# ╟─690df6ec-71ec-46cd-8699-960217bd4f06
# ╟─a0574649-00a6-4883-9a87-0dac3fc5f8a5
# ╟─caaff430-fa09-4bbd-b106-95e36d9203a3
# ╠═90dfe32b-0930-4067-8936-6f1e1e922a35
# ╟─c820cf54-fa0a-4112-8e76-8f76839b7a49
# ╠═673b09d2-5dd4-4b6c-897e-2fc43f04ab8f
# ╟─08b1ac9b-58e8-41de-994f-a05609df3b2c
# ╟─cfb3aaab-2166-43ce-9fdf-586b98fe8272
# ╠═135856f2-6c6b-4bc4-9e7a-ca678d5e729d
# ╠═2f11999f-5d6e-4807-8bce-49791e7a0211
# ╠═ae54aa7a-3afb-43c9-b083-0908b1f02d18
# ╠═d92dc54e-2ba6-4fa1-bd43-7d06c5d9cb6f
# ╠═93a9beb2-084d-4fe7-937c-39d74740cade
# ╠═5941ed8d-62c0-4d11-b897-30fc24f25b78
# ╠═f8f1c3e8-ed99-4a37-8faa-82ec1f946898
# ╠═dde2b892-eae4-4dfd-8735-b533e8a5ae68
# ╠═e4acf92a-3242-4579-a046-95649f835c36
# ╟─835d1acf-95f6-47c1-9bdd-0b0a75034353
# ╠═29497f7b-795e-433f-b772-72191f52dc24
# ╠═2696b1df-e52c-497c-b62f-a0932da6c8a4
# ╠═eb7eb4ca-6088-487c-9017-6b7988188c20
