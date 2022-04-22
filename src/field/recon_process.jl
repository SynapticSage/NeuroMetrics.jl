module recon_process
    
    using ..field
    include("../table.jl")
    using ProgressMeter
    using DataFrames
    using Infiltrator
    using DataStructures

    si = field.operation.selectind
    sk = field.operation.selectkey

    export get_recon_name
    export get_recon_req
    export perform_reconstructions_marginals_and_error
    

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # RUNNING THE PROCESS OF RECONSTRUCTION
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    shortcut_names = OrderedDict("headdir"             => "H",
                                 "currentHeadEgoAngle" => "γ",
                                 "currentPathLength"   => "p",
                                 "stopWell"            => "G")
    𝕀(d) = Dict(zip(values(shortcut_names), keys(shortcut_names)))[d]
    get_recon_name(x, op) = replace(x, "vs("=>"ε(", ","=>") $op ε(")
    function get_recon_req(recon_compare)
        rr = vec([x[i] for x in values(recon_compare), i in 1:2])
        rr= [Set(rr)...]
    end
    """
    Returns the integer dim indices for a prop-string e.g. "x-y"->[1,2]
    """
    function 𝔻(dimstr,dims)
        out = [findfirst(dim.==dims) for dim in split(dimstr,",")]
        if out isa Vector{Nothing}
            @error "Nope! One your dims=$dimstr is wrong. Check for missing syntax (commas)"
        end
        return out
    end
    """
    Returns the remaining dimensions not covered by a prop-string
    """
    𝔻̅(dimstr,dims) = setdiff(1:length(dims), 𝔻(dimstr,dims)) # dims inverse
    """
    Returns remaning dimensions as a joined prop string, instead of ints
    """
    𝔻̅ⱼ(dimstr,dims) = join(dims[𝔻̅(dimstr, dims)],",") # joined
    """
    Returns string with shortcut names
    and the ₀ version returns the remaining names
    """
    ℝ(dims)   = replace(dims, [shortcut_names...][begin:1:end]...) # Replace
    ℝ₀ⱼ(dims,props) = ℝ(join(props[𝔻₀(dims)],"-"))                       # Joined

    function perform_reconstructions_marginals_and_error(beh, spikes,
            K::NamedTuple; recon_compare::Union{Dict, Nothing}=nothing,
            F::Dict=Dict(), recon_summary::DataFrame=DataFrame())

        if recon_compare == nothing
            throw(ArgumentError("Please pass recon_compare"))
        end
        recon_req = get_recon_req(recon_compare)
        @debug "props = $(K.props)"
        local dims  = ℝ(K.props)
        @debug "dims = $(dims)"
        
        x = Set(vec([split(x,"|")[1] for x in recon_req])) # requires the LHS andd inverse of the RHS of each reconstruction
        y = Set(vec([𝔻̅ⱼ(split(x,"|")[2], dims) for x in recon_req])) # requires the LHS andd inverse of the RHS of each reconstruction
        marginals_required = x ∪ y

        @time X = field.get_fields(beh, spikes; K...);
        field_size = size(si(X.Rₕ))

        # ---------
        # MARGINALS
        # ---------
        # Acquire marginals P(X,Y), P(γ, p, G)
        @time @showprogress for marginal in marginals_required
            d̅ =  𝔻̅(marginal, dims)
            println("marginal=>$marginal d̅ = $(d̅)")
            @time F[marginal] = operation.marginalize(X, dims = d̅ );
        end

        # ---------------
        # Reconstructions
        # ---------------
        # Obtain reconstructions!
        R̂ = Dict()
        @time for reconstruction in recon_req
            if reconstruction == nothing
                continue
            end
            @debug "reconstruction = $reconstruction"
            dimr, given = split(reconstruction, "|")
            #inverse_given = join(dims[𝔻₀(given)], ",")
            #ig_set   = split(inverse_given,",")
            @debug "dimr=$dimr, $dims=dims"
            marginalize_dims = 𝔻̅(dimr, dims)
            println(reconstruction)
            @time R̂[reconstruction] = operation.apply(recon.reconstruction, 
                                                      X.occR, 
                                                      F[given].Rₕ;
                                                      marginalize_dims=marginalize_dims,
                                                     );
            @assert ndims(si(R̂[reconstruction]))     == length(split(dimr,","))
            @assert all(size(si(R̂[reconstruction])) .== field_size[𝔻(dimr, dims)])
        end

        # ---------------
        # SUMMARIES
        # ---------------
        # Get reconstruction model error summary
        E = Vector{DataFrame}([])
        for reconstruction in recon_req
            what, given = split(reconstruction, "|")
            error = recon.reconstruction_error(F[what].Rₕsq, R̂[reconstruction])
            error = table.to_dataframe(error; name="error")
            push!(E, error)
        end
        E = vcat(E..., source=:model=>recon_req)
        what,under = [vec(x) for x in eachrow(cat(split.(E.model,"|")...; dims=2))]
        E.what, E.under = what, under
        E = sort(E, [:model, :area, :unit])
        recon_summary = vcat(recon_summary, E)
        field.utils.pushover("Finished reconstruction summaries")
        return recon_summary
    end
end
