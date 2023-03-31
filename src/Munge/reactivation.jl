__precompile__(false)
"""
    reactivation.jl

Contains functions for reactivation analysis.

Detecting cell assemblies in large neuronal populations    
https://www.sciencedirect.com/science/article/pii/S0165027013001489?via%3Dihub
"""
module reactivation

    using MultivariateStats
    using DimensionalData
    using DataFrames
    using StatsBase
    using Infiltrator
    using Missings
    using LinearAlgebra

    export zscoreFRmatrix
    """
        zscoreFRmatrix(X::Union{DimArray,Matrix})

    Returns a matrix of z-scores of the columns of X, which
    are neurons and rows are time bins. 
    """
    function zscoreFRmatrix(X::Union{DimArray,Matrix})
        cols = eachcol(X) 
        cols = map(cols) do x
                    miss = ismissing.(x) .&& .!isnan.(x)
                    x[.!miss] = zscore(disallowmissing(x[.!miss]))
                    x
        end
        hcat(cols...)
    end 
    function zscoreFRmatrix(X::DataFrame)
        cell_columns = names(X)[findall(tryparse.(Int, names(X)) .!== nothing)]
        zscoreFRmatrix(Matrix(X[:, cell_columns]))
    end

    export correlationMatrix
    """
        function correlationMatrix(zX::Union{DimArray,Matrix},
                                   zY::Union{DimArray,Matrix}=zX)::Matrix

    Returns the correlation matrix of z-scores of the firing rates of neurons.
    
    # Arguments
    - `zX::Union{DimArray,Matrix}`: Matrix of z-scores of the firing rates of neurons.
    - `zY::Union{DimArray,Matrix}`: Matrix of z-scores of the firing rates of neurons.

    # Returns
    - `Matrix`: Correlation matrix.
    """
    correlationMatrix(zX::Union{DimArray,Matrix}) = zX' * zX / size(zX, 1)
    correlationMatrix(zX::Union{DimArray,Matrix},
                      zY::Union{DimArray,Matrix}) = zX' * zY / size(zX, 1)

    export getPattern_fromICA
    """
        function getPattern_fromICA(X::Union{DimArray,Matrix};
            maxiter=100, tol=1e-5, do_whiten=true, usefastica=false)::ICA

    Returns the precorrelation matrix from the ICA of X.

    # Arguments
    - `X::Union{DimArray,Matrix}`: Matrix of z-scores of the firing rates of neurons.
    - `maxiter::Int`: Maximum number of iterations.
    - `tol::Float64`: Tolerance.
    - `do_whiten::Bool`: Whether to whiten the data.
    - `usefastica::Bool`: Whether to use fastica.
    
    # Returns
    - `ICA`: ICA object.
    """
    function getPattern_fromICA(C::Union{DimArray,Matrix}, k=nothing;
        maxiter=100, tol=1e-5, do_whiten=true, usefastica=false)
        if k === nothing
            k = size(C, 2)
        end
        ica = fit(ICA, C', k;
            maxiter, tol, do_whiten)
        ica
    end

    """
        estimateKfromPCA(X::Union{DimArray,Matrix})::Int
    Returns the number of components to use for ICA from the PCA of X.
    # Arguments
    - `X::Union{DimArray,Matrix}`: Matrix of z-scores of the firing rates of
                                   neurons.
    """
    function estimateKfromPCA(Z::Union{DimArray,Matrix})
        pca = fit(PCA, Z, maxoutdim=size(Z, 2))
        # Fit a Marcenko-Pastur distribution to the eigenvalues of the covariance
        p_lambda = marcenko_pasteur(eigvals(pca), pca.proj)
        # Find the first eigenvalue that is less than the Marcenko-Pastur distribution
        sum(eigvals(pca) .> p_lambda.lambda_max)
    end

    export marcenko_pasteur
    """
        marcenko_pasteur(lambdas::Vector{Float64}, pca_predict::Matrix)
    Returns the Marcenko-Pasteur distribution for the eigenvalues of the covariance
    # Arguments
    - `lambdas::Vector{Float64}`: Eigenvalues of the covariance matrix.
    - `pca_predict::Matrix`: PCA projection matrix.
    # Returns
    - `Vector{Float64}`: Marcenko-Pasteur distribution.
    https://bnikolic.co.uk/blog/python/2019/11/28/marchenko-pastur.html
    """
    function marcenko_pasteur(lambdas::Vector{Float64}, pca_predict::Matrix)
        @infiltrate
        B = size(pca_predict, 1) # time bins
        N = size(pca_predict, 2) # number of dimensions
        # lamdamax = (1+sqrt((total_ctx_neuron + total_hp_neuron)/length(qW1_spike(1,:)))).^2;
        gamma = B / N # also called q in some papers
        @assert gamma ≥ 1 # has to be greater than one
        sigma² = 1 # because we're using z-scores
        lambda_max = sigma² * (1 + sqrt(gamma))^2
        lambda_min = sigma² * (1 - sqrt(gamma))^2
        m(x) = max.(x, 0)
        p_lambda = ((gamma / (2π * sigma²)) * 
        sqrt.(m(lambda_max .- lambdas) .* m(lambdas .- lambda_min)) ./ lambdas)
        (;lambda_max, lambda_min, p_lambda)
    end

    export reactscore
    """
        function reactivationScore(Zarea1, W, Zarea2)

    Returns the reactivation score between two areas, Zarea1 and Zarea2,

    Zarea1 and Zarea2 are matrices of z-scores of the firing rates of neurons in
    area 1 and area 2, respectively. W is the pattern matrix.

    # Arguments
    - `Zarea1::Matrix`: Matrix of z-scores of the firing rates of neurons in area 1.
    - `W::Matrix`: Project pattern matrix containing the assemblies
    - `Zarea2::Matrix`: Matrix of z-scores of the firing rates of neurons in area 2.
    """
    function reactscore(Zarea1::Matrix, P::Matrix, Zarea2::Matrix)
        r = Zarea1 * P; 
        R = r .* Matrix(Zarea2) 
    end
    """
        function reactivationScore(Zarea1, W, Zarea2)

    Returns the reactivation score between two areas, Zarea1 and Zarea2,
    
    # Arguments
    - `Zarea1::Matrix`: Matrix of z-scores of the firing rates of neurons in area 1.
    - `W::ICA`: ICA object containing the assemblies
    - `Zarea2::Matrix`: Matrix of z-scores of the firing rates of neurons in area 2.
    """
    function reactscore(Zarea1, P::ICA, Zarea2)
        reactscore(Zarea1, P.W, Zarea2)
    end

    """
        reactivationScore(Zarea1, Zarea2, k=nothing)

    Returns the reactivation score between two areas, Zarea1 and Zarea2
    
    # Arguments
    - `Zarea1::Matrix`: Matrix of z-scores of the firing rates of neurons in area 1.
    - `Zarea2::Matrix`: Matrix of z-scores of the firing rates of neurons in area 2.
    - `k::Int`: Number of components to use for ICA.
    """
    function reactscore(Zarea1::Matrix, Zarea2::Matrix, P, 
        ica::Union{Nothing,ICA}; 
            k::Union{Int,Nothing}=nothing)
        if Zarea1 !== Zarea2
            Z, Zarea1, Zarea2 = get_Z(Zarea1, Zarea2)
        end
        if  ica === nothing
            Zarea1 = Matrix(predict(ica, Zarea1')')
            Zarea2 = Matrix(predict(ica, Zarea2')')
            Z = Matrix(predict(ica, Z')')
        end
        reactscore(Zarea1, P, Zarea2)
    end

    """
        function get_Z(Zarea1::Matrix, Zarea2::Matrix)

    Transforms the z-scores of the firing rates of neurons in area 1 and area 2
    into a single matrix, Z, and returns Z, and returns zero-padded single-area
    matrices, Zarea1 and Zarea2.
    """
    function get_Z(Zarea1::Matrix, Zarea2::Matrix)
        if Zarea1 === Zarea2
            Z = Zarea1
        else
            Z = [Zarea1 Zarea2]
            z1 = zeros(size(Zarea1))
            z2 = zeros(size(Zarea2))
            Zarea1 = [Zarea1 z2]
            Zarea2 = [z1 Zarea2]
            @assert size(Zarea1, 2) == size(Zarea2, 2)
        end
        Z = disallowmissing(Z)
        Zarea1 = disallowmissing(Zarea1)
        Zarea2 = disallowmissing(Zarea2)
        Z, Zarea1, Zarea2
    end

    """
        function reactivationSauce(Zarea1::Matrix, Zarea2::Matrix, k=nothing)

    Returns the materials/sauce needed to calculate the reactivation score between
    two areas, Zarea1 and Zarea2, Zarea1 and Zarea2 are matrices of z-scores of the
    firing rates of neurons in area 1 and area 2, respectively.
    """
    function reactivationSauce(Zarea::Matrix, Zarea2::Matrix,
        k=Union{Int,Nothing}=nothing, ica=true)
        Z, Zarea, Zarea2 = get_Z(Zarea, Zarea2)
        if ica
            if k === nothing
                k = estimateKfromPCA(Z)
            end
            ica    = getPattern_fromICA(Z, k)
            P = correlationMatrix(Z)
        else
            ica = nothing
            P = correlationMatrix(Z)
        end
        return (;P, ica)
    end


end
