
# MLJulia based :: fewer cross-validation options

import MLJBase: StratifiedCV, ResamplingStrategy, train_test_pairs
using ScientificTypes
using CategoricalArrays

function get_field_crossval_jl(beh::DataFrame, data::DataFrame, 
        shifts::Union{StepRangeLen,Vector{T}} where T <: Real;
        cv::T=StatifiedCV(nfolds=2) where T <: ResamplingStrategy, 
        stratify_cols::Union{Nothing,Union{String},String}=nothing,
        kws...)

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


function visualize_crossvals(beh::DataFrame, gen)
    tr_ranges = []
    te_ranges = []
    plots = []
    plots1 = []
    for (i,(train, test)) in enumerate(gen)
        push!(tr_ranges, extrema(beh[train,:time]))
        push!(te_ranges, extrema(beh[test,:time]))
        push!(plots,  Plots.histogram(beh.time[train], label="$i train"))
        push!(plots,  Plots.histogram(beh.time[test],  label="$i test"))
        push!(plots1, Plots.plot(beh.time[train],      label="$i train"))
        push!(plots1, Plots.plot(beh.time[test],       label="$i test"))
    end
    Plots.plot(plots...)
    Plots.plot(plots1...)
end




import MLJBase: StratifiedCV, ResamplingStrategy, train_test_pairs
using ScientificTypes
using CategoricalArrays

function get_field_crossval_jl(beh::DataFrame, data::DataFrame, 
        shifts::Union{StepRangeLen,Vector{T}} where T <: Real;
        cv::T=StatifiedCV(nfolds=2) where T <: ResamplingStrategy, 
        stratify_cols::Union{Nothing,Union{String},String}=nothing,
        kws...)

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

# ScikitLearn based

import Sklearn.CrossValidation
function get_field_crossval_sk(beh::DataFrame, data::DataFrame, 
        shifts::Union{StepRangeLen,Vector{T}} where T <: Real;
        cv::T=StatifiedCV(nfolds=2) where T <: ResamplingStrategy, 
        stratify_cols::Union{Nothing,Union{String},String}=nothing,
        kws...)

    equalize = beh[!, stratify_cols]
    equalize = convert.(eltype(equalize), equalize)
    multiple_cols = ndims(equalize) == 1
    equalize = multiple_cols ? Vector(equalize) : Matrix(equalize)
    if multiple_cols
        @error "not implemented yet"
    end

end

