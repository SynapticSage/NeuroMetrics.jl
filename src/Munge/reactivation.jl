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
    - `Z::Union{DimArray,Matrix}`: Matrix of z-scores of the firing rates of neurons.
    - `maxiter::Int`: Maximum number of iterations.
    - `tol::Float64`: Tolerance.
    - `do_whiten::Bool`: Whether to whiten the data.
    - `usefastica::Bool`: Whether to use fastica.
    
    # Returns
    - `ICA`: ICA object.
    """
    function getPattern_fromICA(Z::Union{DimArray,Matrix}, k=nothing;
        maxiter=100, restarts=10, tol=1e-5, do_whiten=true, usefastica=false)
        if k === nothing
            k = estimateKfromPCA(Z)
        end
        # where Z has dimensions (time bins, neurons/dimensions)
        @assert(size(Z, 1) ≥ size(Z, 2), 
                "time bins must be greater than neurons")
        @info "ICA with $k components, maxiter=$maxiter, tol=$tol, do_whiten=$do_whiten" size(Z)
        ica = except = nothing
        while ica === nothing && restarts > 0
            try
                restarts -= 1
                ica=fit(ICA, Z', k; maxiter, tol, do_whiten)
            catch e
                @warn "ICA failed" e restarts
                except = ConvergenceException(maxiter, NaN, tol)
            end
        end
        ica isa ICA ? ica : except
    end

    """
        estimateKfromPCA(X::Union{DimArray,Matrix})::Int
    Returns the number of components to use for ICA from the PCA of X.
    # Arguments
    - `X::Union{DimArray,Matrix}`: Matrix of z-scores of the firing rates of
                                   neurons.
    """
    function estimateKfromPCA(Z::Union{DimArray,Matrix})
        @info "Estimating K from PCA" size(Z)
        pca = nothing
        try
            pca = fit(PCA, Z, maxoutdim=size(Z, 2))
        catch e
            @infiltrate
            throw(e)
        end
        # Fit a Marcenko-Pastur distribution to the eigenvalues of the covariance
        p_lambda = marcenko_pasteur(eigvals(pca), pca.proj)
        # Find the first eigenvalue that is less than the Marcenko-Pastur distribution
        K = sum(eigvals(pca) .> p_lambda.lambda_max)
        @info "K from PCA" K
        K
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

    """
        preICA_Z
    Preprocessed Z-transformed firing rate matrices
    # Fields
    - `Z::Matrix`: Z-transformed firing rate matrix, concatenated if area1 and area2 are different.
    - `Zarea1::Matrix`: Z-transformed firing rate matrix for area 1.
    - `Zarea2::Matrix`: Z-transformed firing rate matrix for area 2.
    """
    struct preICA_Z
        Z::Matrix
        Zarea1::Matrix
        Zarea2::Matrix
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
        Z      = disallowmissing(Z)
        Zarea1 = disallowmissing(Zarea1)
        Zarea2 = disallowmissing(Zarea2)
        # Remove any all NaN columns
        @infiltrate
        Z      = Z[:,findall(.!vec(all(isnan.(Z), dims=1)))]
        Zarea1 = Zarea1[:, findall(.!vec(all(isnan.(Zarea1), dims=1)))]
        Zarea2 = Zarea2[:, findall(.!vec(all(isnan.(Zarea2), dims=1)))]
        # Remove any all NaN rows
        Z      = Z[findall(.!vec(all(isnan.(Z), dims=2))), :]
        Zarea1 = Zarea1[findall(.!vec(all(isnan.(Zarea1), dims=2))), :]
        Zarea2 = Zarea2[findall(.!vec(all(isnan.(Zarea2), dims=2))), :]
        preICA_Z(Z, Zarea1, Zarea2)
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
        if size(Zarea1, 2) == size(P,1)
            r = Zarea1 * P; 
            R = r .* Matrix(Zarea2) 
        else
            r = Zarea1' * P; 
            R = r .* Matrix(Zarea2') 
        end
    end

    export react
    """
        function react(Zarea1::Matrix, Zarea2::Matrix, k=nothing)

    Returns the materials/sauce needed to calculate the reactivation score
    between two areas, Zarea1 and Zarea2, Zarea1 and Zarea2 are matrices of
    z-scores of the firing rates of neurons in area 1 and area 2, respectively.
    # Arguments
    - `Zarea1::Matrix`: Matrix of z-scores of the firing rates of neurons in area 1.
    - `Zarea2::Matrix`: Matrix of z-scores of the firing rates of neurons in area 2.
    - `k::Int`: Number of components to use for ICA.
    - `ica::Bool`: Whether to use ICA or not.
    """
    function react(preZ::preICA_Z;
        k::Union{Int,Nothing}=nothing, ica=true)
        Z, Zarea1, Zarea2 = preZ.Z, preZ.Zarea1, preZ.Zarea2
        if ica
            ica  = getPattern_fromICA(Z, k)
            if typeof(ica) <: ConvergenceException
                P = nothing
            else
                P = correlationMatrix(Matrix(predict(ica, Z')'))
            end
        else
            ica = nothing
            P = correlationMatrix(Z)
        end
        return (;P, ica, k)
    end
    function react(Zarea1::Matrix, Zarea2::Matrix;
        k::Union{Int,Nothing}=nothing, ica=true)
        preZ = get_Z(Zarea1, Zarea2)
        react(preZ, k=k, ica=ica)
    end

    """
        function reactivationScore(Zarea1::Matrix, Zarea2::Matrix, k=nothing)
    Returns the reactivation score between two areas, Zarea1 and Zarea2,
    Zarea1 and Zarea2 are matrices of z-scores of the firing rates of neurons in
    area 1 and area 2, respectively.
    # Arguments
    - `Zarea1::Matrix`: Matrix of z-scores of the firing rates of neurons in area 1.
    - `Zarea2::Matrix`: Matrix of z-scores of the firing rates of neurons in area 2.
    - `k::Int`: Number of components to use for ICA.
    - `ica::Bool`: Whether to use ICA or not.
    """
    function reactscore(Zarea1::Matrix, Zarea2::Matrix;
        k::Union{Int,Nothing}=nothing, ica=true)
        preZ      = get_Z(Zarea1, Zarea2)
        P, ica, k = react(preZ, k=k, ica=ica)
        if typeof(ica) <: ConvergenceException
            @warn "ICA did not converge" ica
            return NaN
        end
        z1, z2 = predict(ica, preZ.Zarea1'), predict(ica, preZ.Zarea2')
        reactscore(z1, P, z2)
    end

end
