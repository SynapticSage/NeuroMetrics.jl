include("./imports_isolated.jl")
include("./pfc_cyclewise.jl")
using GoalFetchAnalysis
using .Munge.isolated


#  GET SMALLER BATCH DATAFRAME FOR TESTING
# bin period 18
df_sub = @subset(df, :bins .>= 17 .&& :bins .<= 19)

#    _  _      ____           _          
# _| || |_   / ___|__ _  ___| |__   ___ 
#|_  ..  _| | |   / _` |/ __| '_ \ / _ \
#|_      _| | |__| (_| | (__| | | |  __/
#  |_||_|    \____\__,_|\___|_| |_|\___|
#                                       
function create_caches_xy(df_sub)
    dx_dy = ThreadSafeDict()
    Threads.@threads for relcyc in collect(-opt["cycles"]:opt["cycles"])
        dx_dy[(;relcyc)] = get_dx_dy(df_sub, relcyc)
    end
    dx_dy = Dict(dx_dy)
end

function create_caches_xyb(df_sub)
    DFB= groupby(df_sub, :bins)
    prog = Progress((opt["cycles"]*2+1) * length(DFB))
    @time dx_dy_bin = ThreadSafeDict()
    sets = [(relcyc, DFB[b], b) 
        for relcyc in -opt["cycles"]:opt["cycles"],
        b in eachindex(DFB)]
    Threads.@threads for (relcyc, dfb, b) in sets
        key = (;relcyc, bin=b[1])
        try
            dx_dy_bin[key] = get_dx_dy(dfb, relcyc)
        catch exception
            @infiltrate
            @warn "$key failed" exception
            sleep(0.1)
        end
        next!(prog)
    end
    dx_dy_bin = Dict(dx_dy_bin)
    return dx_dy_bin
end

dx_dy     = create_caches_xy(df_sub)
dx_dy_bin = create_caches_xyb(df_sub)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sets to iterate GLM over
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = [b for b in collect(keys(dx_dy_bin))]
uAreas = unique(cells.area)
bf(F) = [(b, f) for b in B, f in F]

param_interactions_hc = [ (dist, dep, ind) for dist in ["probit"],
     # (dep, ind) in Iterators.product(("PFC","CA1"),("PFC", "PFC"))]
     (dep, ind) in zip(("CA1",),("PFC",))]
hcF(ind,dep) = construct_predict_isospikecount(df_sub, cells, ind; dep_area=dep)
param_interactions_hc = [ (b, f, dist, dep, ind) 
    for (dist, dep, ind) in param_interactions_hc,
    (b, f) in bf(hcF(ind,dep))]

param_interactions_sc = [ (dist, dep, ind) for dist in ["softplus"],
    # (dep, ind) in Iterators.product(("PFC","CA1"),("PFC", "PFC"))]
     (dep, ind) in zip(("CA1",),("PFC",))]
scF(ind,dep) = construct_predict_spikecount(df_sub, cells, ind; dep_area=dep)
param_interactions_sc = [ (b, f, dist, dep, ind) 
    for (dist, dep, ind) in param_interactions_sc,
    (b, f) in bf(scF(ind,dep))]


begin # test first sets
    (dist, dep, ind) = first(param_interactions_sc)
    (b, f) = first(bf(construct_predict_spikecount(df_sub,    cells, ind; dep_area=dep)))
    (b, f) = first(bf(construct_predict_isospikecount(df_sub, cells, ind; dep_area=dep)))
end
zscoreX = true

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"""
    addconfusion_remove_y_ypred!(shuffle_cellhasiso)
Add a confusion matrix to each model in `shuffle_cellhasiso` and remove `y` and
`ypred` from all but the last N entries.
"""
function addconfusion_remove_y_ypred!(shuffle_cellhasiso;N=100)
    # Add a confusion matrix to each model
    for (k,v) in shuffle_cellhasiso
        if !haskey(v, "y") || !haskey(v, "ypred")
            # @warn "skipping" k typeof(v)
            continue
        end
        y, ypred = v["y"], v["ypred"]
        cm = MLJ.ConfusionMatrix()(y, ypred)
        v["cm"] = cm
    end
    # Except for the last N entries shuffle, remove y and ypred from the model
    # dicts
    for (k,v) in Iterators.take(shuffle_cellhasiso, 
                                length(shuffle_cellhasiso)-N)
        if !haskey(v, "y") || !haskey(v, "ypred")
            continue
        end
        delete!(v, "y")
        delete!(v, "ypred")
        if haskey(v, "yr")
            delete!(v, "yr")
            delete!(v, "ypredr")
        end
    end
    # Confusion matrix of the last shuffle y, ypred and yr, ypredr
    shuf = collect(values(shuffle_cellhasiso))[end]
    y, ypred = shuf["y"], shuf["ypred"]
    yr, ypredr = shuf["yr"], shuf["ypredr"]
    cm = MLJ.ConfusionMatrix()(y, ypred)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_spikecount   = initorget("model_spikecount"; obj=Dict())
shuffle_spikecount = initorget("shuffle_spikecount"; obj=Dict())
summary(model_spikecount); summary(shuffle_spikecount); 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    """
        postprocess_model!(dx_dy_bin, dx_dy_key, f::FormulaTerm, dist; 
                shuffle=false,
                model::Dict, other_keys...)

    Postprocess the model `model` for the given `dx_dy_key` and `f` and `dist`.
    # Arguments
    - `dx_dy_bin`: the `dx_dy_bin` dictionary
    - `dx_dy_key`: the key for the `dx_dy_bin` dictionary
    - `f`: the `FormulaTerm` used to create the model
    - `dist`: the distribution used to create the model
    - `shuffle`: whether the model was created by shuffling the data
    
    # Return
    - `key`: the key for the `model` dictionary
    - `model[key]`: the model
    """
    function poisson_model!(dx_dy_bin, dx_dy_key, f::FormulaTerm, dist; 
                shuffle=false,
                model::Dict, other_keys...)
        XXb, yb = dx_dy_bin[dx_dy_key]
        if shuffle
            # under shuffle, XXb and yb might be of different sizes
            m=min(size(XXb,1), size(yb,1))
            XXb, yb = XXb[1:m, :], yb[1:m, :]
        end
        unit=parse(Int, replace(string(f.lhs)))
        key = (;other_keys..., dx_dy_key..., unit)
        if !opt["overwrite"] && haskey(model, key)
            return key, nothing
        end
        nonmiss = vec(any( (!).(ismissing.(Array(yb)) .|| 
                    ismissing.(Array(XXb))), dims=2))
        XXb, yb = XXb[nonmiss, cellcols(XXb)], 
                   yb[nonmiss, cellcols(yb)]
        for i in axes(XXb)[2]
            XXb[!,i] = replace(XXb[:,i], missing=>0)
        end
        for i in axes(yb)[2]
            yb[!,i] = replace(yb[:,i], missing=>0)
            # yb[:, i] = Int.(yb[!, i] .> 0)
        end
        if ismissing(yb) || ismissing(XXb)
            return key, nothing
        end
        mod = glm_(f, XXb, yb, dist; type=:pyglm,  zscoreX=zscoreX,
            cv=2, tol=1e-2, max_iter=50,
            desc=Dict("bin"=>key.bin, "unit"=>key.unit, "relcyc"=>key.relcyc,
                      "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
                      "dep"=>dep, "ind"=>ind))
        model[key] = mod
        if haskey(mod, "score") && !isnan(mod["score"])
            # @exfiltrate; 
            # @warn "key=$key has score=$(mod["score"])" mod["y"], mod["ypred"]
        end
        return key, mod
    end

    # =============================================
    # PREDICT SPIKCOUTNS ALL TIMES BY 10 minute BIN
    #
    # Distribution: POISSON / SOFTPLUS
    # =============================================
    @showprogress "single run" for (b, f, dist, dep, ind) in param_interactions_sc
        XXb, yb = dx_dy_bin[b]
        answer = poisson_model!(dx_dy_bin, b, f, dist; model=model_spikecount, 
                                        dir="$ind => $dep", dist)
        if answer[2] === nothing
            @warn "skipping" answer[1]
            continue
        end
        # @infiltrate
    end

    commit_cycwise_vars("model_spikecount")

    # ------------ Shuffle ------------------------
    opt["shuffle"] = 100

    # dx_dy_bin = deepcopy(dx_dy_bin)
    @showprogress "SHUFFLE" for _ in 1:opt["shuffle"]
        shuf = hash(rand())
        println("SHUFFLE=$shuf shuffling...")
        dx_dy_bin_shuffle = create_caches_xyb(isolated.shuffle_cyclelabels(df_sub))
        println("SHUFFLE=$shuf shuffling...done")
        @showprogress "single run" for (b, f, dist, dep, ind) in 
        param_interactions_sc
            answer = poisson_model!(dx_dy_bin_shuffle, b, f, dist; 
                    shuffle=true,
                    model=shuffle_spikecount, shuf, dir="$ind => $dep", dist)
            if answer[2] === nothing
                @warn "skipping" answer[1]
                continue
            end
        end
    end

    commit_cycwise_vars("shuffle_spikecount")
    # ------------ Shuffle ------------------------


    #    _  _     _   _             _           
    #  _| || |_  | | | | __ _ ___  (_)___  ___  
    # |_  ..  _| | |_| |/ _` / __| | / __|/ _ \ 
    # |_      _| |  _  | (_| \__ \ | \__ \ (_) |
    #   |_||_|   |_| |_|\__,_|___/ |_|___/\___/ 
    #                                           
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    """
    logistic_model!(dx_dy_key, f::FormulaTerm; 
                model::Dict, other_keys...)
    Creates a single logistic regression model for a given
    `dx_dy_key` and `f` and stores it in the `model` dictionary.
    # Arguments
    - `dx_dy_key` is a NamedTuple of the form
            (bin=1, relcyc=1)
            which indexes into the Main scope's `dx_dy_bin` 
            a dictionary that caches XX, y to run GLM on
    - `f` is a FormulaTerm that is used to construct the GLM
    - `model` is a dictionary that is used to cache the results
    - `other_keys` are additional keys that are used to index into the 
            `model` dictionary
    # Returns
    - `key` is a NamedTuple that is used to index into the `model` dictionary
    - `m` is the GLM model dict item
    """
    function logistic_model!(dx_dy_bin, dx_dy_key, f::FormulaTerm, dist; 
                shuffle=false,
                model::Dict, other_keys...)
        XXb, yb = dx_dy_bin[dx_dy_key]
        if shuffle
            # under shuffle, XXb and yb might be of different sizes
            m=min(size(XXb,1), size(yb,1))
            XXb, yb = XXb[1:m, :], yb[1:m, :]
        end
        unit=parse(Int, replace(string(f.lhs), "_i"=>""))
        key = (;other_keys..., dx_dy_key..., unit)
        if !opt["overwrite"] && haskey(model, key)
            return key, nothing
        end
        nonmiss = vec(any( (!).(ismissing.(Array(yb)) .|| 
                    ismissing.(Array(XXb))), dims=2))
        XXb, yb = XXb[nonmiss, cellcols(XXb)], yb[nonmiss, 
                               isolatedcellcols(yb)]
        for i in axes(XXb)[2]
            XXb[!,i] = replace(XXb[:,i], missing=>0)
        end
        for i in axes(yb)[2]
            yb[!,i] = replace(yb[:,i], missing=>0)
            yb[:, i] = Int.(yb[!, i] .> 0)
        end
        if ismissing(yb) || ismissing(XXb)
            return key, nothing
        end
        mod = glm_(f, XXb, yb, dist; type=:pyglm,  zscoreX=zscoreX,
            cv=2, tol=1e-2, max_iter=50,
            desc=Dict("bin"=>key.bin, "unit"=>key.unit, "relcyc"=>key.relcyc,
                      "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
                      "dep"=>dep, "ind"=>ind))
        model[key] = mod
        if haskey(mod, "score") && !isnan(mod["score"])
            # @exfiltrate; 
            # @warn "key=$key has score=$(mod["score"])" mod["y"], mod["ypred"]
        end
        return key, mod
    end

    # =============================================
    # PREDICT CELL HAS ISO SPIKE BY 10 minute BIN
    #
    # Distribution: BINOMIAL / PROBIT
    # =============================================
    model_cellhasiso = Dict() #initorget("model_cellhasiso", Dict())
    @showprogress "single run" for (b, f, dist, dep, ind) in param_interactions_hc
        answer = logistic_model!(dx_dy_bin, b, f, dist; model=model_cellhasiso, 
                                        dir="$ind => $dep", dist)
        if answer[2] === nothing
            @warn "skipping" answer[1]
            continue
        end
        # @infiltrate
    end

    commit_cycwise_vars("model_cellhasiso")

    # Filter out the values that have NaN score
    model_cellhasiso = Dict(k=>v for (k,v) in model_cellhasiso if 
        haskey(v, "score") && !isnan(v["score"]))


    # ------------ Shuffle ------------------------
    shuffle_cellhasiso = Dict() #initorget("shuffle_spikecount", Dict())

    opt["shuffle"] = 100
    # dx_dy_bin = deepcopy(dx_dy_bin)
    @showprogress "SHUFFLE" for _ in 1:opt["shuffle"]
        shuf = hash(rand())
        println("SHUFFLE=$shuf shuffling...")
        dx_dy_bin_shuffle = create_caches_xyb(isolated.shuffle_cyclelabels(df_sub))
        println("SHUFFLE=$shuf shuffling...done")
        @showprogress "single run" for (b, f, dist, dep, ind) in 
                 param_interactions_hc
            answer = logistic_model!(dx_dy_bin_shuffle, b, f, dist; 
                    shuffle=true,
                    model=shuffle_cellhasiso, shuf, dir="$ind => $dep", dist)
            if answer[2] === nothing
                @warn "skipping" answer[1]
                continue
            end
            # @infiltrate
        end
    end


    addconfusion_remove_y_ypred!(shuffle_cellhasiso)
    commit_cycwise_vars("shuffle_cellhasiso")
    # shuffle_cellhasiso = initorget("shuffle_cellhasiso"; obj=Dict())

