module recon_process
    
    using ..Field
    import Table
    using ProgressMeter
    using DataFrames
    using Infiltrator
    using DataStructures
    using ThreadSafeDicts
    using NaNStatistics

    si = Field.operation.selectind
    sk = Field.operation.selectkey

    export get_recon_name
    export get_recon_req
    export perform_reconstructions_marginals_and_error
    export var_shortcut_names
    

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # RUNNING THE PROCESS OF RECONSTRUCTION
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    """
    shortcut_names

    shortcut single char names for common behavioral variables
    """
    shortcut_names = OrderedDict("headdir"             => "H",
                                 "currentHeadEgoAngle" => "γ",
                                 "currentPathLength"   => "p",
                                 "stopWell"            => "G")
    const var_shortcut_names = shortcut_names


    get_shortcutnames(items)  = [replace(item, recon_process.var_shortcut_names...)
                             for item in items]
    """
    UnTranslate from shorcut names
    """
    inv_shortcutnames(items)  = [replace(item, Dict(kv[2]=>kv[1] for
                                 kv in shortcut_names)...) for item in items]

    """
    𝕀

    invert a dictionary mapping
    """
    𝕀(d) = Dict(zip(values(shortcut_names), keys(shortcut_names)))[d]

    """
    get_recon_name

    get the name we prefer to use for reconstruction errors
    """
    get_recon_name(x, op) = replace(x, "vs("=>"ε(", ","=>") $op ε(")

    """
    get_recon_req

    return all of the required reconstructions for a set of comparisons
    """
    function get_recon_req(comparisons::AbstractDict)::Vector{String}
        rr = vec([x[i] for x in values(comparisons), i in 1:2])
        rr= [Set(rr)...]
    end

    """
    𝔻 : string2num

    Returns the integer dim indices for a prop-string e.g. "x-y"->[1,2]
    """
    function 𝔻(dimstr,dims)
        out = [findfirst(dim.==dims) for dim in split(dimstr,",")]
        if out isa Vector{Nothing}
            @error "Nope! One your dims=$dimstr is wrong. Check for missing syntax (commas)"
        end
        return out
    end
    const string2num = 𝔻

    """
    𝔻̅ : string2numinv

    Returns the remaining dimensions not covered by a prop-string
    """
    𝔻̅(dimstr,dims) = setdiff(1:length(dims), 𝔻(dimstr,dims)) # dims inverse
    const string2numinv = 𝔻̅

    """
    𝔻̅ⱼ : string2stringinv
    Returns remaning dimensions as a joined prop string, instead of ints
    """
    𝔻̅ⱼ(dimstr,dims) = join(dims[𝔻̅(dimstr, dims)],",") # joined
    const string2stringinv = 𝔻̅ⱼ

    """
    Returns string with shortcut names
    and the ₀ version returns the remaining names
    """
    ℝ(dims)   = replace(dims, [shortcut_names...][begin:end]...) # Replace
    ℝ₀ⱼ(dims,props) = ℝ(join(props[𝔻₀(dims)],"-"))                 # Joined

    function perform_reconstructions_marginals_and_error(beh, spikes,
            K::NamedTuple; recon_compare::Union{Dict, Nothing}=nothing,
            R̂::AbstractDict=ThreadSafeDict(),
            F::AbstractDict=Dict(), recon_summary::DataFrame=DataFrame())

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

        @time X = Field.get_fields(beh, spikes; K...);
        field_size = size(si(X.Rₕ))

        # ---------
        # MARGINALS
        # ---------
        # Acquire marginals P(X,Y), P(γ, p, G)
        P = Progress(length(marginals_required), desc="Marginals")
        @time for marginal in marginals_required
            d̅ =  𝔻̅(marginal, dims)
            println("marginal=>$marginal d̅ = $(d̅)")
            @time F[marginal] = operation.marginalize(X, dims = d̅ );
            next!(P)
        end

        # ---------------
        # Reconstructions
        # ---------------
        # Obtain reconstructions!
        P = Progress(length(recon_req), desc="Reconstruction")
        @time @Threads.threads for reconstruction in recon_req
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
            next!(P)
        end

        # ---------------
        # SUMMARIES
        # ---------------
        # Get reconstruction model error summary
        E = Vector{DataFrame}([])
        P = Progress(length(recon_req), desc="Reconstruction error")
        for reconstruction in recon_req
            what, given = split(reconstruction, "|")
            error = recon.reconstruction_error(F[what].Rₕsq, R̂[reconstruction])
            error = Table.to_dataframe(error; name="error")
            push!(E, error)
            next!(P)
        end
        E = vcat(E..., source=:model=>recon_req)
        what,under = [vec(x) for x in eachrow(cat(split.(E.model,"|")...; dims=2))]
        E.what, E.under = what, under
        E = sort(E, [:model, :area, :unit])
        recon_summary = vcat(recon_summary, E)
        Field.utils.pushover("Finished reconstruction summaries")
        return recon_summary
    end

    # ----------------
    # Munging
    # ----------------
    function create_unstacked_error_table(E, recon_compare)
        uE = unstack(E[!,Not([:what,:under])], :model,:error)[!,Not([:dim_1, :dim_2])]
        uE.∑ε = vec(nansum(replace(Matrix(uE[:, get_recon_req(recon_compare)]),missing=>NaN); dims=2))
        uE = sort(uE, [:area,:∑ε])
        for (rc, compare) in recon_compare
            #uE[!,rc*"div"] = uE[!, compare[1]] ./ uE[!, compare[2]]
            println(rc)
            uE[!,rc] = (uE[!, compare[1]] .- uE[!, compare[2]])./
                       (uE[!,compare[1]] .+ uE[!,compare[2]])
            uE[!,rc] = (uE[!, compare[1]] .- uE[!, compare[2]])
        end
        uE
    end
end
