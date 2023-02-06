module crossval

    export crossval
    export get_field_crossval_sk, get_field_crossval_jl

    # MLJulia based :: fewer cross-validation options
    import MLJBase: StratifiedCV, ResamplingStrategy, train_test_pairs
    using ScientificTypes
    using CategoricalArrays
    using DataFrames
    import DIutils
    using Table: DFColVars

    in_range = DIutils.in_range


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

    function perform_split(train::Vector, test::Vector, beh::DataFrame, data::DataFrame)

        beh.ind  = 1:size(beh,1)
        data.ind = 1:size(data,1)

        beh.type .= ""
        beh.type[train, :] .= "train"
        beh.type[test, :]  .= "test"

        beh, data = raw.register(beh, data; on="time", transfer="type")
        nonmatchingepochs = (!).(DIutils.ismember(data.epoch, unique(beh.epoch)))
        data.type[nonmatchingepochs] .= ""

        return beh[beh.type .== "train", :], 
               beh[beh.type .== "test", :],
               data[data.type .== "train", :], 
               data[data.type .== "test", :]
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
            shuffle::Bool=true,
            stratify_cols::Union{Nothing,Union{String},String}=nothing,
            return_cv::Bool=true,
            kws...)

        if cv !== CrossValidation.KFold
            G = DIutils.findgroups((beh[:,col] for col in stratify_cols)...)
            if cv === CrossValidation.StratifiedKFold
                @info "Stratified K-fold random=$shuffle"
                cv = cv(G; n_folds, shuffle)
            elseif cv === CrossValidation.LabelKFold
                @info "Labeled K-fold random=$shuffle"
                cv = cv(G; n_folds, shuffle)
            end
        else
            @info "K-fold random=$shuffle"
            cv = cv(size(beh,1); n_folds, shuffle)
        end

        results = Dict()
        for (f,(train, test)) in enumerate(cv)
            B_train, B_test, D_train, D_test = perform_split(train, test,
                                                             beh, data)
            test_key, train_key = _tk("test",  f), _tk("train", f)
            results[test_key] = timeshift.get_field_shift(B_train, D_train, shifts;
                                                          kws...)
            results[train_key] = timeshift.get_field_shift(B_test, D_test, shifts;
                                                           kws...)
        end

        if return_cv
            return results, cv
        else
            return results
        end
    end

    function unpack_cv_as_df(cvresult; shift_scale=:minutes)
        K = keys(cvresult)
        df_k = DIutils.namedtupkeys_to_df(K)
        n_folds = maximum(df_k.fold)
        result = DataFrame()
        for fold in 1:n_folds
            train = timeshift.info_to_dataframe(cvresult[_tk("train", fold)], shift_scale=shift_scale)
            test  = timeshift.info_to_dataframe(cvresult[_tk("test", fold)], shift_scale=shift_scale)
            train.fold .= fold
            test.fold .= fold
            combined = vcat(train, test; 
                           cols=:intersect, 
                           source=:source=>[:train, :test])
            append!(result, combined)
        end
        result
    end

    using Metrics: r2_score, binary_accuracy, categorical_accuracy

    """
    since we don't have ground truth, we compare the distribution
    or point estimate found in train-indices to those in test-indices
    """

    function apply_metric_to_traintest(cvresult::AbstractDict;
            groupby_vars::DFColVars,
            stackid_vars::DFColVars,
            metric=r2_score, shift_scale=:minutes)

            combined = unpack_cv_as_df(cvresult; shift_scale)
            apply_metric_to_traintest(combined; groupby_vars, stackid_vars, 
                                      metric)
    end

    function apply_metric_to_traintest(combined::DataFrame;
            groupby_vars::DFColVars,
            stackid_vars::DFColVars,
            metric=mae)

        combined = dropmissing(unstack(combined, [groupby_vars,stackid_vars],
                                      :source, :info))

        combined = groupby(combined, groupby_vars)

        for C in combined
            key = first(C[!,groupby_vars])
            if occursin("adjusted", String(Symbol(metric)))
                M = metric(C.train, C.test, size(C,1))
            else
                M = metric(C.train, C.test)
            end
            if key isa DataFrameRow
                D = DataFrame(key)
                D[!, :fold] = fold
                D[!, Symbol(metric)] = M
            else
                D = DataFrame(groupby_vars=>key, 
                              Symbol(metric)=>M,
                              :fold=>fold)
            end
            push!(result, D)
        end

        result = vcat(result...)

        return result
    end

    function visualize_crossval()
    end

end
