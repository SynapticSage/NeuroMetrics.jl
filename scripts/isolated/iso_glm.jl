include("./imports_isolated.jl")
include("./pfc_cyclewise.jl")

#  GET SMALLER BATCH DATAFRAME FOR TESTING
# bin period 18
df_sub = @subset(df, :bins .>= 17 .&& :bins .<= 19)

#    _  _      ____           _          
# _| || |_   / ___|__ _  ___| |__   ___ 
#|_  ..  _| | |   / _` |/ __| '_ \ / _ \
#|_      _| | |__| (_| | (__| | | |  __/
#  |_||_|    \____\__,_|\___|_| |_|\___|
#                                       
function create_caches(df_sub)
    dx_dy = ThreadSafeDict()
    Threads.@threads for relcyc in collect(-opt["cycles"]:opt["cycles"])
        dx_dy[(;relcyc)] = get_dx_dy(df_sub, relcyc)
    end
    dx_dy = Dict(dx_dy)
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
    return dx_dy, dx_dy_bin
end
dx_dy, dx_dy_bin = create_caches(df_sub)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sets to iterate GLM over
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = [b for b in collect(keys(dx_dy_bin))]
uAreas = unique(cells.area)
bf(F) = [(b, f) for b in B, f in F]
param_interactions_sc = [ (dist, dep, ind) for dist in ["poisson", "softplus"],
    (dep, ind) in (("PFC","CA1"),("PFC", "PFC"))]
param_interactions_hc = [ (dist, dep, ind) for dist in ["binomial", "probit"],
     (dep, ind) in (("PFC","CA1"),("PFC", "PFC"))]
begin # test first sets
    (dist, dep, ind) = first(param_interactions_sc)
    (b, f) = first(bf(construct_predict_spikecount(df_sub,    cells, ind; dep_area=dep)))
    (b, f) = first(bf(construct_predict_isospikecount(df_sub, cells, ind; dep_area=dep)))
end
zscoreX = true

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_spikecount   = initorget("model_spikecount"; obj=Dict())
shuffle_spikecount = initorget("shuffle_spikecount"; obj=Dict())
summary(model_spikecount); summary(shuffle_spikecount); 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    # =============================================
    # PREDICT SPIKCOUTNS ALL TIMES BY 10 minute BIN
    #
    # Distribution: POISSON / SOFTPLUS
    #
    # =============================================
    @showprogress "param-area-interact" for (dist, dep, ind) in 
                param_interactions_sc
        F = construct_predict_spikecount(df_sub, cells, ind; dep_area=dep)
    @showprogress "single run" for (b,f) in bf(F)
        XXb, yb = dx_dy_bin[b]
        for i in 1:size(XXb,2)
            XXb[!,i] = replace(XXb[:,i], missing=>0)
        end
        for i in 1:size(yb,2)
            yb[!,i] = replace(yb[:,i], missing=>0)
        end
        key = (;b..., unit=parse(Int,string(f.lhs)), dir="$ind => $dep", 
                dist)
        if !opt["overwrite"] && haskey(model_spikecount, key)
            continue
        end
        m = glm_(f, XXb, yb, dist; type=:pyglm,  zscoreX=zscoreX,
            cv=2, tol=1e-2, max_iter=50,
            desc=Dict("bin"=>key.bin, "unit"=>key.unit, "relcyc"=>key.relcyc,
                      "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
                    "dist"=>dist, "dep"=>dep, "ind"=>ind))
        model_spikecount[key] = m
    end
    end
    commit_cycwise_vars("model_spikecount")

    # ------------ Shuffle ------------------------
    # shuffle_spikecount = initorget("shuffle_spikecount")
    # @showprogress "shuffle" for _ in 1:opt["shuffle"]
    #     shuf = hash(rand())
    # @showprogres "area-interact" for (dep, ind) in ind_dep
    #     F = construct_predict_spikecount(df_sub, cells, ind; dep_area=dep)
    # @showprogress "ones shuf" for (b,f) in bf(F)
    #     XXb, yb = dx_dy_bin[b]
    #     key = (;shuf, b..., unit=parse(Int,string(f.lhs)))
    #     if !opt["overwrite"] && haskey(shuffle_spikecount, key)
    #         continue
    #     end
    #     s = shuffle_glm_(f, XXb, yb, "poisson"; type=:pyglm,
    #         desc=Dict("bin"=>key.bin, "unit"=>key.unit, "relcyc"=>key.relcyc,
    #                   "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
    #             "dep"=>dep, "ind"=>ind))
    #     shuffle_spikecount[key] = s
    # end
    # end
    # end
    # commit_cycwise_vars("shuffle_spikecount")
    # ------------ Shuffle ------------------------



    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model_cellhasiso   = initorget("model_cellhasiso"; obj=Dict())
    shuffle_cellhasiso = initorget("shuffle_cellhasiso"; obj=Dict())
    summary(model_cellhasiso); summary(shuffle_cellhasiso)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # =============================================
    # PREDICT CELL HAS ISO SPIKE BY 10 minute BIN
    #
    # Distribution: BINOMIAL / PROBIT
    #   
    # =============================================
    @showprogress "param-area-interact" for (dist, dep, ind) in 
                param_interactions_hc
        F = construct_predict_isospikecount(df_sub, cells, ind; dep_area=dep)
    @showprogress "single run" for (b,f) in bf(F)
        XXb, yb = dx_dy_bin[b]
        nonmiss = vec(any(
                 (!).(ismissing.(Array(yb)) .|| ismissing.(Array(XXb)) ),
            dims=2))
        XXb, yb = XXb[nonmiss, cellcols(XXb)], yb[nonmiss, isolatedcellcols(yb)]
        for i in 1:size(XXb,2)
            XXb[!,i] = replace(XXb[:,i], missing=>0)
        end
        for i in 1:size(yb,2)
            yb[!,i] = replace(yb[:,i], missing=>0)
            yb[:, i] = Int.(yb[!, i] .> 0)
        end
        unit=parse(Int, replace(string(f.lhs), "_i"=>""))
        b = (;b..., unit, dir="$ind => $dep", dist)
        if ismissing(yb) || ismissing(XXb) ||
            (!opt["overwrite"] && haskey(model_cellhasiso, b))
            continue
        end
        m = glm_(f, XXb, yb, dist; type=:pyglm,  zscoreX=zscoreX,
            cv=2, tol=1e-2, max_iter=50,
            desc=Dict("bin"=>b.bin, "unit"=>b.unit, "relcyc"=>b.relcyc,
                      "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
                "dep"=>dep, "ind"=>ind))
        model_cellhasiso[b] = m
        sleep(0.001)
    end
    end
    commit_cycwise_vars("model_cellhasiso")

    # ------------ Shuffle ------------------------
    # shuffle_cellhasiso = initorget("shuffle_spikecount")
    # @showprogrss "hyperparmas" for dist in ["binomial", "probit"]
    # @showprogress "shuffle" for _ in 1:opt["shuffle"]
    #     shuf = hash(rand())
    # @showprogres "area-interact" for (dep, ind) in ind_dep
    #     F = construct_predict_spikecount(df_sub, cells, ind; dep_area=dep)
    # @showprogress "ones shuf" for (b,f) in bf(F)
    #     XXb, yb = dx_dy_bin[b]
    #     key = (;shuf, b..., unit=parse(Int,string(f.lhs)))
    #     if !opt["overwrite"] && haskey(shuffle_spikecount, key)
    #         continue
    #     end
    #     s = shuffle_glm_(f, XXb, yb, "poisson"; type=:pyglm,
    #         desc=Dict("bin"=>key.bin, "unit"=>key.unit, "relcyc"=>key.relcyc,
    #                   "unit"=>string(f.lhs), "dir"=>"$ind => $dep",
    #                   "dep"=>dep, "ind"=>ind))
    #     shuffle_spikecount[key] = s
    # end
    # end
    # end
    # end
    # commit_cycwise_vars("shuffle_cellhasiso")
    # ------------ Shuffle ------------------------

# =========================
#  . .     ,---.|    ,-.-.
# -+-+-    |  _.|    | | |
# -+-+-    |   ||    | | |
#  ` `     `---'`---'` ' '
#  OF HAS isolated spike per cell  🔺
# ========================
# model_cellhasiso, cacheiso, shuffle_cellhasiso = 
#                             initorget("model_cellhasiso"), ThreadSafeDict(),
#                             initorget("shuffle_cellhasiso"; obj=Dict())
# model_cellhasiso = OrderedDict()
# pos = (df, dx_dy, cells, construct_predict_isospikecount, :pyglm)
# pos = (df, Dict(first(dx_dy)), cells, construct_predict_isospikecount, :pyglm)
# kws = (Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""],
#     ytrans=x->Float64(x>0), xtrans=x->Float64(x), modelz=model_cellhasiso,
#     )
# isolated.run_glm!(pos...;kws...) 

# Relative cycle 0 runs of the model
# isolated.run_glm!(df, dx_dy, cells, construct_predict_spikecount, :pyglm;
#     Dist=Distributions.Binomial(), unitwise=true, unitrep=["_i"=>""],
#     ytrans=x->Float64(x>0), xtrans=x->Float64(x), modelz=model_cellhasiso,
#     )


# # Running for entire df (not by relcyc)
# isolated.run_glm!(df, cells, construct_predict_spikecount, :pyglm;
#     Dist=Distributions.Binomial(), unitwise=true, modelz=model)


# tmp = shuffle_cellhasiso
# shuffle_cellhasiso = ThreadSafeDict()
# for (k,v) in tmp
#     push!(shuffle_cellhasiso,k=>v)
# end


# -----------------------------------------------
# Description of the different construct_ methods
# -----------------------------------------------
#
# construct_predict_isospikecount :: 
#     Predicts the number of isolated spikes in a cycle
# Type :: Poisson
#
# construct_predict_iso :: 
#     Whether there is an isolated spike in a cycle
# Type :: Binomial
#
# construct_predict_spikecount :: 
#     Predicts the number of spikes in a cycle
# Type :: Poisson
#
# -----------------------------------------------
# Description of different model variables stored
# in the local JLD2 file
# -----------------------------------------------
#
# model_cellhasiso :: 
#     Model for predicting whether there is an isolated spike in a cycle
#
# model_cellhasiso_matlab ::
#     Model for predicting whether there is an isolated spike in a cycle
#     using the matlab implementation
#
# 
# shuffle_cellhasiso :: 
#     Shuffled model for predicting whether there is an isolated spike in a cycle
# 
# model_spikecount :: 
#     Model for predicting the number of spikes in a cycle
#
# model_counthasiso :: 
#     Model for predicting whether there is an isolated spike in a cycle
#     using the number of spikes in a cycle as a predictor
# 
#

# run_shuffle!(pos...; kws..., 
#     shuffle_models=shuffle_cellhasiso,
#     shufcount=30)

commit_cycwise_vars()
