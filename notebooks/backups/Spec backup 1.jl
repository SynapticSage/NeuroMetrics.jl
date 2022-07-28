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

# ╔═╡ dfcfa875-d122-49f1-ab24-66c1937b3134
begin
	using Revise
	import DrWatson
	DrWatson.quickactivate(expanduser("~/Projects/goal-code"))
	push!(LOAD_PATH, DrWatson.srcdir())
end

# ╔═╡ a9b4b3d2-f318-11ec-210a-a70a7964ee72
Field, Load, Timeshift, Utils, Load, Table = begin
	using GoalFetchAnalysis
	import Timeshift
	import Field
	import Utils
	import Load
	import Table
	Field, Load, Timeshift, Utils, Load, Table
end

# ╔═╡ 450738b2-3d49-4e45-9a4d-ffa1721f833a
begin
	using StatsBase
	using Term
	using REPL.TerminalMenus
	using Plots
	import Arrow
	using DataFramesMeta
	using DataFrames
	using ColorSchemes
	using StatsBase
	using Infiltrator
	using Dates
	using StatsPlots
	using Logging

	if Utils.in_range(hour(now()), [0,5]) ||
	   Utils.in_range(hour(now()), [20, 24])
		Plots.theme(:dark)
		theme="dark"
	else
		Plots.theme(:bright)
		theme="bright"
	end
	
	"Importing packages, theme=$theme"
end

# ╔═╡ 435f9680-1520-468e-b97c-2ea4fb2c1ff4
using PlutoUI

# ╔═╡ 8d41c178-16ee-4881-b55c-fb80f274d7bc
PlutoUI.TableOfContents(title="Non-adaptive field shifting")

# ╔═╡ 42ea762b-12ed-4eb8-ade0-3bffff593690
md"""
# Package Imports 📦
"""

# ╔═╡ 10a1552f-ffda-4b79-8a4c-fe2864bc3ae5
md"""
Import some key GoalFetchAnalysis codes
"""

# ╔═╡ 12c97814-1d83-4f82-9f5f-891abb878e60
import Plot

# ╔═╡ bf904f0e-7387-4b73-890b-fbb5fc6137de
md"""
# Loading data  💾
First, we're going to loadup a key dataframe of interest: **cells**
"""

# ╔═╡ cadaf555-3b90-4c5b-846b-686ce4130497
cells = Load.load_cells("RY16", 36, "*")


# ╔═╡ Cell order:
# ╠═8d41c178-16ee-4881-b55c-fb80f274d7bc
# ╟─dfcfa875-d122-49f1-ab24-66c1937b3134
# ╟─42ea762b-12ed-4eb8-ade0-3bffff593690
# ╟─450738b2-3d49-4e45-9a4d-ffa1721f833a
# ╟─10a1552f-ffda-4b79-8a4c-fe2864bc3ae5
# ╠═a9b4b3d2-f318-11ec-210a-a70a7964ee72
# ╟─12c97814-1d83-4f82-9f5f-891abb878e60
# ╟─bf904f0e-7387-4b73-890b-fbb5fc6137de
# ╟─cadaf555-3b90-4c5b-846b-686ce4130497
# ╠═13b0d82f-d721-40f7-b2f6-d099d41e8897
# ╟─77039e5f-1c59-463d-8299-13dd6af9631c
# ╟─7ee46c7e-eeb4-4ae2-bc4b-318322aff9fb
# ╟─3c02c040-2fed-4bed-b985-8e78bb241455
# ╟─eb62a3a9-5ba8-4737-8214-06b3c436eac2
# ╟─edcb03f0-3e52-4fd6-adbe-48d37150ba13
# ╟─61d024da-011a-42a3-b456-19475da19e78
# ╟─5bad660a-016b-4e05-83c3-38e10e3323a1
# ╟─435f9680-1520-468e-b97c-2ea4fb2c1ff4
# ╠═b553e927-d6b8-469a-90de-b1b0bf9efa11
# ╠═225323c9-4ed6-42ce-987d-4d5557efaa35
# ╟─f77875f8-21c4-4d97-8103-cdc7d33adee3
# ╟─6399b851-1f0a-4999-ae00-e66437f5e264
# ╟─d47c3d1d-394b-45bb-9d8b-1dc8ddc1b9b3
# ╠═bfcf1a0a-6520-4f79-b6c2-f97d8d0f34cc
# ╠═ee772af0-a435-4390-9fe5-d4cd3a2aacac
# ╟─3c5fad0a-fd7f-4766-88fc-53c2ad7bcca4
# ╟─d1ae7695-1e2f-4d90-9663-37f500bfd53a
# ╟─07daf23c-24dc-4907-bdad-ba2c91978f8e
# ╟─35ea2991-0a66-47a6-bec2-d1f44725bc8e
# ╟─131fb039-7631-429b-b327-73a40e408b59
# ╟─5adb69a3-ab98-401e-a344-38aada960e6d
# ╟─3cb16c9d-eb21-4947-8597-3991917cc7f0
# ╟─0713b24c-cdc4-435f-a8c1-872a071b7c50
# ╟─75a5d90c-8b21-4bd9-abd2-1c2794fe31e3
# ╟─84a0f7ed-a1f2-4bf7-82a1-150be1acba53
# ╟─b4f43b70-1ff9-49fd-ad05-1209e4713c39
# ╟─0a9df174-91a1-45fe-93c2-cb787ef6c4e9
# ╟─ec76028d-61eb-447e-bfc6-7a2d46abd411
# ╟─72fdc63f-ce5d-42a7-8d51-b3bb4cbd7624
# ╟─4539b9f8-4ba1-4576-82ca-02bc2e8151d1
# ╟─5ffb3ee8-141d-4ee0-8021-021c746b8bc5
# ╟─4f71df4b-e203-403e-bf52-611f58f9751c
# ╟─3f946154-9c2a-4a31-bba3-99f3e68e0dea
# ╟─44f6cfdd-fd5d-4757-b497-e6b58039d95c
# ╟─6bcc8054-4f6d-4d82-a34a-b5541235d87a
# ╟─33c64062-002a-4eea-8e29-e8e36784a666
# ╟─c8800b8a-fd8b-43b2-a1fa-bbb76879e56d
# ╟─2a75cc20-9bf4-4d2e-9cd5-a5d30d8d8c76
# ╟─e81ceff1-a29f-428c-b448-d3f1ac204ef4
# ╟─63237edf-5fb0-44f7-8650-134974d71734
# ╟─7bcaad0e-17ab-4908-b3f2-0d56b03b0c87
# ╟─d84ddfac-bb8b-4203-a734-7ee0500201b6
# ╟─a557bb0a-a88a-4519-a235-ae446630bdd7
# ╟─765e0db8-1364-437e-8ecf-3d5055d8a07f
# ╟─17b37173-4fbe-49e9-8d80-503ab8d73f97
# ╟─6f425636-86d6-4314-a0f6-5f5772999918
# ╟─410c533c-744a-4ac6-8d6e-1e5a55b2ed27
# ╟─3e0ae61d-c7e5-477f-9576-b76462b7085b
# ╟─25e863d7-6530-4e9d-a482-a06b3379e53f
# ╠═f8ef3758-951d-4d64-96c2-7856caf21fce
# ╟─f46b7018-e10d-4d96-a572-85fd8df371dd
# ╠═13bda800-9e07-4322-a283-fda5aa2869eb
# ╠═0729d86b-3dae-4e40-a2fd-a078c3665277
# ╟─975a6d70-ae6b-4b25-8fef-34552c39c94c
# ╟─b2ed7800-80f1-4d30-a00a-c76bcd8036a9
# ╟─07b01f6f-f0b8-4f23-ac64-53f7dc9ea6f2
# ╟─1762e086-7f10-47cd-a895-fae9a772d6d5
# ╟─09819531-982a-4ded-9b0c-5187aac26e97
# ╟─063b8c00-e52b-4888-a384-b15df02afa8d
# ╟─d728c9d6-4571-40bb-b1df-b34d7ae0b785
# ╟─b5a4f445-1cd2-476a-aeab-a362a2de3b17
# ╟─32ddb00f-1eac-491f-a321-3d9f48c7f70c
# ╟─125653c1-55de-4410-ad30-d00433008005
# ╟─07d3adea-b413-461d-a68d-0383cd5ab26b
# ╟─8ffb6013-7a53-4ff7-804f-89a71cbc7397
# ╟─f2c7d4d1-d673-4d17-af1c-339373dba7c2
# ╟─21080e4d-614d-4b70-8c8d-e70120c9fb02
# ╟─735ab0d6-0c9f-499d-b62d-93b7fcc3b39e
# ╟─3b6aaa97-aee1-46cf-b4cf-79d7312056eb
# ╟─c42b9721-f30f-45b0-a6ab-60413eee9876
# ╠═cd3925de-dbb0-4e57-9a1e-48bf8fbb109f
# ╠═2bb18fdc-9080-4f0a-9d09-c2cbe0e6404a
# ╠═398a837c-4710-4eb6-9f78-786d7173bc49
# ╠═49bb383c-b3cf-448f-9906-ab1eabae3f75
# ╠═54d550c9-0e83-4cdb-bb27-88d4c6dffe83
# ╠═5df33da9-bb8d-4356-99b0-a4742a20c87e
