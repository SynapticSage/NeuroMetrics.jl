
# MLJulia based :: fewer cross-validation options

import MLJBase: StratifiedCV, ResamplingStrategy, train_test_pairs
using ScientificTypes
using CategoricalArrays

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
        cv = cv(1:size(beh,1); n_folds)
    end

    for (train, test) in cv
        B_train, B_test = beh[train,:], beh[test,:]
        strain = utils.in_range(spikes.time, extrema(B_train.time))
        stest  = utils.in_range(spikes.time, extrema(B_test.time))
        D_train, D_test = data[strain,:], data[stest,:]
        results[(;type = "train", fold=f)] = timeshift.get_field_shift(B_train,
                                                                       D_train,
                                                                       shifts; 
                                                                       kws...)
        results[(;type = "test",  fold=f)] = timeshift.get_field_shift(B_test, 
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

function visualize_crossval()
end

