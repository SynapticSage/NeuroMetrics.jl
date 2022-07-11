### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ fc614ab8-00cb-11ed-0f62-f751ef056b39
begin
	
	using DrWatson, Revise
	quickactivate(expanduser("~/Projects/goal-code"))
	using DataFrames
	using KernelDensity, Distributions
	using Plots, Measures
	using ProgressMeter, ProgressLogging
	using StatsPlots
	using DataFramesMeta
	using DataStructures: OrderedDict
	using Distributions
	import Base.Threads: @spawn
	using ThreadSafeDicts, NaNStatistics
	using Combinatorics: powerset
	

	using GoalFetchAnalysis
	import Timeshift
	using Timeshift.dataframe: info_to_dataframe
	import Load
	import Filt
	filts = Filt.get_filters()


    function keymessage(I::AbstractDict, key)
        docontinue=false
        if key ∈ keys(I)
            if I[key] isa Task && !(istaskfailed(I[key]))
                @info "task key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                docontinue=true
            elseif I[key] isa Task && istaskfailed(I[key])
                "key=$key already exists, but failed...redo!"
            else
                @info "key=$key already exists"
                printstyled("SKIPPING...\n", blink=true)
                docontinue=true
            end
        end
        if key ∉ keys(I)
            @info "key=$key ∉ keys, ...creating..."
        end
		docontinue
     end
     
end

# ╔═╡ 823b1bff-d922-4c2b-8a50-179af24094bd
begin
	@time spikes, beh, ripples, cells = Load.load("RY16", 36);
	_, spikes = Load.register(beh, spikes; transfer=["velVec"], on="time")
		
end

# ╔═╡ 90dfe32b-0930-4067-8936-6f1e1e922a35
begin
    if isfile(Timeshift.mainspath())
        I = Timeshift.load_mains()
    else
        I = OrderedDict()
    end
end

# ╔═╡ 673b09d2-5dd4-4b6c-897e-2fc43f04ab8f
begin

    @progress "Datacut iteration" for datacut ∈ collect(keys(filts))
        finished_batch = false
        @progress "Props" for props ∈ prop_set
    #        marginal = 𝕄(props)
    #        key = get_key(;marginal, datacut, shifts)
    #        if keymessage(I, key); continue; end
    #        I[key] = @time Timeshift.get_field_shift(beh, spikes, shifts; newkws...)
    #        finished_batch = true
    #    end
    #    if finished_batch
    #        Timeshift.save_mains(I)
        end
    end

end

# ╔═╡ 135856f2-6c6b-4bc4-9e7a-ca678d5e729d


# ╔═╡ dde2b892-eae4-4dfd-8735-b533e8a5ae68


# ╔═╡ df63982b-553c-4043-b29b-320b440b9883


# ╔═╡ 835d1acf-95f6-47c1-9bdd-0b0a75034353


# ╔═╡ Cell order:
# ╠═fc614ab8-00cb-11ed-0f62-f751ef056b39
# ╠═823b1bff-d922-4c2b-8a50-179af24094bd
# ╠═90dfe32b-0930-4067-8936-6f1e1e922a35
# ╠═673b09d2-5dd4-4b6c-897e-2fc43f04ab8f
# ╠═135856f2-6c6b-4bc4-9e7a-ca678d5e729d
# ╠═dde2b892-eae4-4dfd-8735-b533e8a5ae68
# ╠═df63982b-553c-4043-b29b-320b440b9883
# ╠═835d1acf-95f6-47c1-9bdd-0b0a75034353
