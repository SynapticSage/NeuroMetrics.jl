module timeshift

    __revise_mode__ = :evalassign
    using DataFrames
    import ..field
    export get_field_shift
    export field
    using ThreadSafeDicts
    using DataStructures
    using DataFramesMeta
    using ProgressMeter
    using Distributions
    using Revise
    include("../table.jl")
    include("../shuffle.jl")
    using Infiltrator
    using Distributed
    using Dagger
    using Serialization
    using DrWatson
    using StatsPlots: @df
    using Plots
    using NaNStatistics

    export get_field_shift_shuffles, get_field_shift
    export to_dataframe, info_to_dataframe
    export fetch_best_fields
    export fetch_best_fields
    export shuffle_correct
    export ts_plotdir
    export plot_shifts
    export info_dataframe_and_cell_dataframe
    export save_mains, save_shuffles

    # -------------------- SHIFTING TYPES ---------------------------
    shift_func(data::DataFrame, shift::Real) = 
             transform(data, :time => (t->t.+shift) =>:time, copycols=false)
    const σ = shift_func

    # -------------------- SHIFTED Receptive Fields --------------------------
    function get_field_shift(beh::DataFrame, data::DataFrame, shift::Real; kws...)
        print("Using single")
        fieldobj = field.get_fields(σ(beh, shift), data; kws...)
    end

    function get_field_shift(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            multithread::Bool=true,
            postfunc::Union{Function,Nothing}=nothing,
            as::Type=OrderedDict,
            multi::Union{Bool, Symbol}=true,
            safe_dict::AbstractDict=ThreadSafeDict(),
            kws...)

        kws = (;dokde=false, kws...)

        if multi isa Bool
            multi = multi ? :thread : :single
        end
        @info "Starting multi=$multi"
        msg = "$multi timeshift-shuffles"

        p = Progress(length(shifts), desc="field shift calculations")
        if multi == :thread
            Threads.@threads for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = field.get_fields(σ(beh, shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        elseif multi == :single
            @showprogress for shift in shifts
                if shift ∈ keys(safe_dict)
                    continue
                end
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, shift=>result)
                next!(p)
            end
        elseif multi == :distributed
            @error "Not implemented"
        end
        safe_dict = Dict(safe_dict...)
        out = as(key=>pop!(safe_dict, key) for key in sort([keys(safe_dict)...]))
        return out
    end

    function get_field_shift_shuffles(beh::DataFrame, data::DataFrame,
            shifts::Union{StepRangeLen,Vector{T}} where T <: Real; 
            multi::Union{Symbol,Bool}=true,
            nShuffle::Int=100, 
            shuffle_pos::Union{Tuple,NamedTuple}=(;),
            shuffle_kws::NamedTuple=(;),
            shuffle_func=shuffle.by,
            postfunc::Union{Function,Nothing}=nothing,
            safe_dict::AbstractDict=ThreadSafeDict(),
            exfiltrateAfter::Real=Inf,
            kws...)::AbstractDict

        kws = (;dokde=false, kws...)

        sets = collect( Iterators.enumerate(Iterators.product(1:nShuffle,
                                                              shifts)))

        # Generate the distribution functoin for all shuffles (expensive to do for each)
        if isempty(shuffle_pos) || !(shuffle_pos[1] isa Distribution)
            shuffle_kws = (; shuffledist_df=beh, shuffle_kws...)
            distribution = shuffle._create_distribution(data, 
                                                        shuffle_pos...;
                                                        shuffle_kws...)
            shuffle_pos = (;shuffle_pos..., distribution)
            @info "distribution=$distribution"
        end
        
        if multi isa Bool
            multi = multi ? :thread : :single
        end

        @info "Starting multi=$multi"
        msg = "$multi timeshift-shuffles"
        skips = 0
        P = Progress(length(sets), desc=msg)
        if multi == :thread
            Threads.@threads for (i, (shuffle,shift)) in sets
                if (;shift,shuffle) ∈ keys(safe_dict)
                    skips+=1
                    continue
                end
                data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, (;shift,shuffle)=>result)
                next!(P)
                if mod(i-skips, exfiltrateAfter) == 0
                    @info "chechpoint->exfiltrated"
                    @exfiltrate
                end
            end
        elseif multi == :distributed
            @assert nprocs() > 1 msg
            @showprogress 0.1 msg for (i,(shuffle,shift)) in sets
                data = Dagger.spawn(shuffle_func, data, shuffle_pos...; shuffle_kws...)
                @debug "Dagger 1"
                result = Dagger.spawn(field.get_fields, σ(beh,shift), data; kws...)
                @debug "Dagger 2"
                if postfunc != nothing
                    result = Dagger.@spawn postfunc(result)
                    @debug "Dagger 3"
                end
                push!(safe_dict, (;shift,shuffle)=>result)
            end
        elseif multi == :single
            for (i,(shuffle,shift)) in sets
                if (;shift,shuffle) ∈ keys(safe_dict)
                    skips+=1
                    next!(P)
                    continue
                else
                    #@info (; i,shift,shuffle)
                end
                data = shuffle_func(data, shuffle_pos...; shuffle_kws...)
                result = field.get_fields(σ(beh,shift), data; kws...)
                if postfunc != nothing
                    result = postfunc(result)
                end
                push!(safe_dict, (;shift,shuffle)=>result)
                if mod(i-skips, exfiltrateAfter) == 0
                    @info "chechpoint->exfiltrated"
                    @exfiltrate
                end
                next!(P)
            end
        else
            throw(ArgumentError("Unrecognized argument multi=$multi"))
        end
        safe_dict = Dict(safe_dict...)
        out = OrderedDict(key=>pop!(safe_dict, key) 
                          for key in sort([keys(safe_dict)...]))
        return out
    end


    function to_dataframe(shifts::AbstractDict{<:Real, <:Any}; kws...) 
        table.to_dataframe(shifts; key_name="shift", kws...)
    end
    function to_dataframe(shifts::AbstractDict{<:NamedTuple, <:Any};
            kws...)
        @debug "to_dataframe"
        table.to_dataframe(shifts; kws...)
    end

    function info_to_dataframe(shifts::AbstractDict{<:Real,<:Any};
            kws...)::DataFrame
        table.to_dataframe(shifts; key_name="shift", name="info", kws...)
    end
    function info_to_dataframe(
            shifts::AbstractDict{<:Union{NamedTuple,AbstractDict},<:Any};
            kws...)::DataFrame
        table.to_dataframe(shifts; name="info", kws...)
    end
    function info_dataframe_and_cell_dataframe(place; save_cell_table="", shift_scale=:seconds, kws...)
        df = table.to_dataframe(place, key_name="shift", name="info", kws...)
        df.shift = df.shift .* -1; # beh.time-shift... negative shift is future, so correcting this
        df = if shift_scale == :minutes
            @info "minutes"
            df = transform(df, :shift => (x-> x*60) => :shift)
        elseif shift_scale != :seconds
            @error "Must be seconds or minutes"
        else
            df
        end
        df   = sort(df, [:area,:unit,:shift])
        taus = unique(sort(df.shift))
        df_imax = combine(groupby(df, [:unit, :area]), 
                          :info=>argmax, 
                          :info=>(x->taus[argmax(x)])=>:bestTau)
        df_imax = df_imax[df_imax.bestTau.!=taus[1],:] # Excluding samples with the first tau, because that's the null condition for no variation
        return df, df_imax
    end


    function fetch_best_fields(fieldInfo::DataFrame, pos...; kws...)
        beh, data = pos
        X = Dict()
        for neuron in fieldInfo.units
            D = @subset(data, :unit.==neuron)
            x = get_fields(σ(beh,fieldInfo.bestTau), D; kws...)
            push!(X,x)
        end
        return X
    end

    function func(main, shuffle, transform::Union{Symbol,Pair}; 
            outer=[], by=[:shuffle,:unit], metric=:info, stat=:median, op=.-)
            
        stat      isa Symbol ? stat = eval(stat) : stat
        transform isa Symbol : transform => transform : transform

        shuffle = combine(groubpy(shuffle, by), metric => stat => main)
        by = [b for b in by if b in names(main)]
        by = [by..., outer...]
        main, shuffle = groupby(main, by), groupby(shuffle, by)
        result = DataFrame()
        inputfields, outputfields = transform
        for (m, s) in zip(main, shuffle)
            r = m[!, Not(statopfields)]
            r[!,outputfields] = op(m[!,], stat(s[!,statopfields]))
            append!(result, r)
        end
    end
    function significant(main, shuffle, transform::Union{Symbol, Pair}; kws...)
        kws = (;kws..., op = ((a,b) -> mean(a .> b , dims=1)))
        func(main, shuffle, transform; kws...)
    end
    function correction(main, shuffle, transform; stat=:mean, kws...)
        kws = (;kws..., op=.-, stat=stat)
        func(main, shuffle, transform; kws...)
    end

    # ------------------
    # SAVING AND LOADING
    # ------------------

    """
    _pathshiftdat

    provides path for saved data

    TODO make more flexible
    """
    function _pathshiftdat(;metric=nothing, shifts=nothing, 
            tag="", props=nothing, splitby=nothing, kws...)

        if props isa Nothing
            @error "field keyword props must be given"
        end
        if splitby isa Nothing
            @error "field keyword splitby must be given"
        end

        parent_folder = datadir("exp_pro", "timeshift")
        if !(isdir(parent_folder))
            mkdir(parent_folder)
        end

        tag = isempty(tag) ? tag : "_$tag"
        if shifts != nothing
            start, stop = round(shifts[begin],digits=3),
                          round(shifts[end],  digits=3)
            N = length(shifts)
            shifts = "_Nstartstop=($N,$start:$stop)"
        else
            @error "No shifts name provided"
        end
        if metric == nothing
            @warn "No metric name provided, assuming metric='field'"
            metric = "field"
        end
        jf(x) = join(x,'-')
        props   = "props=$(jf(props))"
        splitby = "_splitby=$(jf(splitby))"
        name = joinpath(parent_folder, "$props$splitby$shifts$tag.serial")
    end

    function ts_plotdir(;file="", kws...)
        parent_folder = plotsdir("timeshift")
        if !(isdir(parent_folder))
            mkdir(parent_folder)
        end
        name = replace(basename(_pathshiftdat(;kws...)), ".serial"=>"")
        name = joinpath(parent_folder, name)
        if !(isdir(name))
            mkdir(name)
        end
        joinpath(name,file)
    end

    function saveshifts(main=nothing, shuffle=nothing; overwrite=false, kws...)

        name = _pathshiftdat(;kws...)
        if isfile(name)
            D = deserialize(name)
        else
            D = Dict()
        end

        d = Dict(:main     => main,
                 :shuffle  => shuffle,
                 :shifts   => kws.shifts,
                 :metric   => kws.metric,
                 :fieldkws => fieldkws)
        D = overwrite ? merge(D, d) : merge(d, D)
        serialize(name, D)
    end

    function save_mains(M::AbstractDict)
        parent_folder = datadir("exp_pro", "timeshift")
        name = joinpath(parent_folder, "mains")
        if isfile(name)
            D = deserialize(name)
        else
            D = Dict()
        end
        D = merge(D, M)

        serialize(name, D)
    end
    function save_shuffles(S::AbstractDict)
        parent_folder = datadir("exp_pro", "timeshift")
        serialize(joinpath(parent_folder, "mains"), S)
    end

    function loadshifts(;kws...)::Dict
        name = _pathshiftdat(;kws...)
        D = deserialize(name)
    end

    # --------
    # PLOTTING
    # --------

    function plot_shifts(place; desc="", shift_scale=:minutes, clim=:cell)

        # List of desirables
        # ------------------
        # - shuffle plots
        #   - all of the below plots with shuffle
        #       - correction
        #       - probability
        # - unify time annotations across plots
        # - normalized version of plots :: where information is [0:lowest for neuron,
        #                                                       1:highest for neuron]
        # - quantile based clims
        #   currently clipping may be changing my gut feeling about the data
        # - save the imax dataframe

        descSave = replace(desc, ":"=>"", " "=>"-")
        function saveplotshift(args...)
            savefig(plotsdir(args[1:end-1]..., args[end] * ".png"))
            savefig(plotsdir(args[1:end-1]..., args[end] * ".svg"))
            savefig(plotsdir(args[1:end-1]..., args[end] * ".pdf"))
        end

        # DENSITY ALL CELLS
        df, df_imax = info_dataframe_and_cell_dataframe(place; shift_scale)
        xlabel = "shift\n(seconds)"
        xlim = (nanminimum(df.shift), nanmaximum(df.shift))

        @df df_imax density(:bestTau, group=:area, 
                     title="$desc BestTau(Information)", xlim=xlim,
                     xlabel=xlabel, ylabel="Density")
        saveplotshift("fields", "shifts",
             "$(descSave)_density_x=shift,y=densBestTau_by=area")
        @df df_imax histogram(:bestTau; group=:area,  xlim,
                     title="$desc BestTau(Information)", 
                     xlabel=xlabel, ylabel="Density")
        saveplotshift("fields", "shifts",
             "$(descSave)_histogram_x=shift,y=densBestTau_by=area")

        # HISTOGRAM ALL CELLS
        df_m = combine(groupby(df, [:area, :shift]), :info=>mean)
        @df df_m bar(:shift, :info_mean; group=:area, xlim,
                     title="$desc Median(Information)", 
                     xlabel=xlabel, ylabel="Shannon MI Bits")
        saveplotshift("fields", "shifts",
             "$(descSave)_histogram_x=shift,y=medianinfo_by=area.svg")

        # PER CELL
        df = sort(df, [:shift,:area,:unit])
        df_u = sort(unstack(df, :shift, :info), [:area, :unit])
        shifts = parse.(Float32,names(df_u)[3:end])
        units = df_u.unit
        areas = df_u.area
        area_divide = findfirst([diff(areas .== "PFC"); 0].==1)
        exclude = Not([x for x in [:area,:unit,:taus]
                       if String(x) in names(df_u)])
        df_u.taus = vec([x[2] for x in argmax(Matrix(df_u[!,exclude]),dims=2)])
        df_u = sort(df_u, [:area, :taus])
        function get_area(area) 
            # cells x shifts
            M = Matrix(df_u[df_u.area.==area, Not([:area,:unit, :taus])])
            M = replace(M, missing=>NaN)
            M[findall(vec([x[1] for x in any(M.!=0,dims=2)])),:]
        end
        norm₁(x, m, M) = begin
            #@infiltrate nanmaximum(x) > 0
            (max.(min.(x,M),m) .- m)./(M-m)
        end
        ca1, pfc = get_area("CA1"), get_area("PFC")
        if clim == :area
            ca1_clim = nanquantile.(vec(ca1), [0.01, 0.95])
            pfc_clim = nanquantile.(vec(pfc), [0.01, 0.95])
        elseif clim == :cell
            ca1_clim = hcat(nanquantile(ca1, 0.01, dims=2),
                            nanquantile(ca1, 0.95, dims=2))
            pfc_clim = hcat(nanquantile(pfc, 0.01, dims=2),
                            nanquantile(pfc, 0.95, dims=2))
            for (i,cc) in zip(1:size(ca1,1), eachrow(ca1_clim))
                ca1[i,:] = norm₁(ca1[i,:], cc[1], cc[2])
            end
            for (i, pp) in zip(1:size(pfc,1), eachrow(pfc_clim))
                pfc[i,:] = norm₁(pfc[i,:], pp[1], pp[2])
            end
            ca1_clim = (0,1)
            pfc_clim = (0,1)
        else
        end
        p = Plots.plot(
                   heatmap(shifts, 1:size(get_area("CA1"),1), ca1, 
                    clims=ca1_clim, colorbar_title="MI (bits)", colorbar_titlefontrotation=0),
                   heatmap(shifts, 1:size(get_area("PFC"),1), pfc, 
                    clims=pfc_clim, colorbar_title="MI (bits)", colorbar_titlefontrotation=0),
                   title="$desc\nMI", xlabel=xlabel, ylabel="cell"
        )
        vline!(p[1], [0], c=:white, linestyle=:dash, label="Zero lag")
        vline!(p[2], [0], c=:white, linestyle=:dash, legendposition=:none)
        saveplotshift("fields", "shifts", "$(descSave)_heatmap_x=shift,y=cellsort_by=area.pdf")
    end

end
