
# MLJulia based :: fewer cross-validation options

import MLJBase: StratifiedCV, ResamplingStrategy, train_test_pairs
using ScientificTypes
using CategoricalArrays

in_range = utils.in_range

function _tk(type::String, fold::Int)
    (;type, fold)
end

function get_field_crossval_jl(beh::DataFrame, data::DataFrame, 
        shifts::Union{StepRangeLen,Vector{<:Real}};
        cv::T=StatifiedCV(nfolds=2), 
        stratify_cols::Union{Nothing,Union{String},String}=nothing,
        kws...) where T <: ResamplingStrategy

    equalize = beh[!, stratify_cols]
    equalize = convert.(eltype(equalize), equalize)
    multiple_cols = ndims(equalize) == 1
    equalize = multiple_cols ? Vector(equalize) : Matrix(equalize)
    if multiple_cols
        @error "not implemented yet"
    end
    equalize = categorical(equalize)

    gen = train_test_pairs(cv, 1:size(beh,1), equalize)
end

function get_field_crossval_jl(beh::DataFrame, data::DataFrame, 
        shifts::Union{StepRangeLen,Vector{<:Real}};
        cv::T=StatifiedCV(nfolds=2), 
        stratify_cols::Union{Nothing,Union{String},String}=nothing,
        kws...) where T <: ResamplingStrategy

    equalize = beh[!, stratify_cols]
    equalize = convert.(eltype(equalize), equalize)
    multiple_cols = ndims(equalize) == 1
    equalize = multiple_cols ? Vector(equalize) : Matrix(equalize)
    if multiple_cols
        @error "not implemented yet"
    end
    equalize = categorical(equalize)
    gen = train_test_pairs(cv, 1:size(beh,1, kws...), equalize)
end

# ScikitLearn based
# ==================
# * Options
# ** KFold : K contiguous chunks
# ** LabelKFold: Pick chunks of labeled rows with matching labels. Rows with
# matching labels move together, as if they are the atomic units (keep labels
# together!)
# ** StratifiedFold: Balance samples of the group/label/y, so that each fold
# contains an equal fraction of each label (split labels apart!)

import ScikitLearn: CrossValidation
using .CrossValidation: KFold
function get_field_crossval_sk(beh::DataFrame, data::DataFrame, 
        shifts::Union{StepRangeLen,Vector{<:Real}};
        n_folds::Int=2, cv::Function=KFold, 
        stratify_cols::Union{Nothing,Union{String},String}=nothing,
        return_cv::Bool=true,
        kws...)

    if cv !== CrossValidation.KFold
        G = utils.findgroups((beh[:,col] for col in stratify_cols)...)
        if cv === CrossValidation.StratifiedKFold
            @info "Stratified K-fold"
            cv = cv(G; n_folds)
        elseif cv === CrossValidation.LabelKFold
            @info "Labeled K-fold"
            cv = cv(G; n_folds)
        end
    else
        @info "K-fold"
        cv = cv(size(beh,1); n_folds)
    end

    results = Dict()
    for (f,(train, test)) in enumerate(cv)
        B_train, B_test = beh[train,:], beh[test,:]
        strain = in_range(data.time, extrema(B_train.time))
        stest  = in_range(data.time, extrema(B_test.time))
        D_train, D_test = data[strain,:], data[stest,:]
        results[_tf("train", f)] = timeshift.get_field_shift(B_train,
                                                                       D_train,
                                                                       shifts; 
                                                                       kws...)
        results[_tk("test",  f)] = timeshift.get_field_shift(B_test, 
                                                                       D_test,
                                                                       shifts;
                                                                       kws...)
    end

    if return_cv
        return results, cv
    else
        return results
    end
end

using Metrics
"""
since we don't have ground truth, we compare the distribution
or point estimate found in train-indices to those in test-indices
"""
function apply_metric_to_traintest(cvresult, metric=mae)
    K = keys(cvresult)
    df_k = utils.namedtupkeys_to_df(K)
    n_folds = maximum(df_k.fold)
    result = []
    for fold in n_folds
        train = cvresult[_tk("train", f)]
        test  = cvresult[_tk("test", f)]
        R = field.operation.binary(train, test, op=metric)
        push!(result, R)
    end
end

function visualize_crossval()
end

