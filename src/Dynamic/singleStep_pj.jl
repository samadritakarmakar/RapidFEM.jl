#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#




#=This algorithm is as per the algorithm mentioned in The Finite Element
Method: Its Basis and Fundamentals, Sixth edition, O.C. Zienkiewicz
Ch-17 The time dimension – discrete approximation in time=#

"""As per the mentioned reference, the A matrix is calcualated as,

A = θ[p-2]*Δt^(p-2)/factorial(p-2)*M + θ[p-1]*Δt^(p-1)/factorial(p-1)*K + θ[p]*Δt^(p)/factorial(p)*C

There is a difference between the reference formula in Zienkiewicz's book and the one used since a θ[0] = 1 is not
possible in julia and θ[1] = 1 is taken rest of elements are the same as θ_Array.
"""
function build_A(M_C_K_Array::Array{SparseMatrixCSC{Float64,Int64},1}, Δt::Float64, thetaModified::Array{Float64, 1})
    @assert length(M_C_K_Array) <= length(thetaModified) "Length of θ_array must be less 1 than, equal to, or greater than M_C_K_Array"
    A::SparseMatrixCSC = spzeros(size(M_C_K_Array[1])...)

    p::Int64 = length(thetaModified) - length(M_C_K_Array) + 1
    for n ∈ 1:length(M_C_K_Array)
        A += (thetaModified[p]*Δt^(p-1)/float(factorial(p-1)))*M_C_K_Array[n]
        p += 1
    end
    return A
end

"""This function evaluates the values of ̇ḟ, f̈, ..., etc. For this it uses taylor series and the values of fₙ, fₙ₊₁, fₙ₊₂, etc."""
function taylor_f_Dot_eval(f_Array::Array{Array{Float64,1},1}, Δt::Float64)
    rowlength::Int64 = length(f_Array)-1
    f_n_Matrix::Array{Float64, 2} = Array{Float64, 2}(undef, rowlength, length(f_Array[1]))
    for i ∈ 1:rowlength
        f_n_Matrix[i,:] = (f_Array[i+1]-f_Array[1])'
    end
    Δt_Matrix::Array{Float64, 2} = Array{Float64, 2}(undef, rowlength, rowlength)
    Δt_Matrix[:,1] = Δt .*collect(1:rowlength)
    for i ∈ 2:rowlength
        Δt_Matrix[:,i] = (Δt_Matrix[:,1].^i)/float(factorial(i))
    end
    f_DotMatrix::Array{Float64, 2} = Δt_Matrix\f_n_Matrix
    #println("f_DotMatrix = ",f_DotMatrix)
    return f_DotMatrix
end

"""This function evaluates the mean values of ̇the solution array, u_mean, u̇_mean, ü_mean, etc """
function buildSolMeanArray(SolutionArray::Array{Array{Float64, 1},1}, Δt::Float64,
    θ_Array::Array{Float64,1})

    SolMeanArray::Array{Array{Float64, 1},1} = Array{Array{Float64, 1},1}(undef, length(SolutionArray))
    for dot ∈ 1:length(SolutionArray)
        p::Int64 = 1
        SolMeanArray[dot] = zeros(length(SolutionArray[1]))
        SolMeanArray[dot] = SolutionArray[dot]
        for i ∈ dot+1:length(SolutionArray)
            #println(p)
            SolMeanArray[dot] += θ_Array[p]*Δt^(p)/factorial(p)*SolutionArray[i]
            p +=1
        end
        #p = 0
    end
    #println("SolMeanArray = ", SolMeanArray)
    return SolMeanArray
end

"""Builds f_mean so the solution can be evaluated as α = A⁻¹f_mean"""
function build_f_mean(f_Array::Array{Array{Float64,1},1}, Δt::Float64, θ_Array::Array{Float64, 1})
    @assert length(f_Array) == length(θ_Array)+1 "Length of f_Array must be 1 greater than length of θ_Array"
    f_DotMatrix::Array{Float64, 2} = taylor_f_Dot_eval(f_Array, Δt)
    f_mean::Array{Float64, 1} = deepcopy(f_Array[1])
    factorial_i::Array{Float64, 1} = Array{Float64, 1}(undef, size(f_DotMatrix,1))
    for i ∈ 1:size(f_DotMatrix,1)
        factorial_i[i] = float(factorial(i))
    end
    for j ∈ 1:size(f_DotMatrix,2)
        for i ∈ 1:size(f_DotMatrix,1)
            f_mean[j] += θ_Array[i]*Δt^(i)/factorial_i[i]*f_DotMatrix[i,j]
        end
    end
    #println("f_mean = ", f_mean)
    return f_mean
end

"""Return A matrix and f_mean vector so that the form 𝐀α + 𝐌ü_mean + 𝐂u̇_mean + 𝐊u_mean + 𝐟_mean = 0 can be formed.

    A::SparseMatrixCSC, SolMeanArray::Array{Array{Float64, 1},1}, f_mean::Array{Float64, 1} = get_SSpj_A_u_f(M_C_K_Array, SolutionArray, f_Array, Δt, θ_Array)
"""
function get_SSpj_A_meanU_f(M_C_K_Array::Array{SparseMatrixCSC{Float64,Int64},1}, SolutionArray::Array{Array{Float64, 1},1},
    f_Array::Array{Array{Float64,1},1}, Δt::Float64, θ_Array::Array{Float64,1})

    thetaModified::Array{Float64, 1} = Array{Float64, 1}(undef, length(θ_Array)+1)
    thetaModified[1] = 1.0
    thetaModified[2:end] = θ_Array
    A::SparseMatrixCSC = build_A(M_C_K_Array, Δt, thetaModified)
    f_mean::Array{Float64, 1} = build_f_mean(f_Array, Δt, θ_Array)
    SolMeanArray::Array{Array{Float64, 1},1} = buildSolMeanArray(SolutionArray, Δt, θ_Array)
    return A, SolMeanArray, f_mean
end

"""Returns Sparse Matrix A and vector b so that α = A⁻¹b can be evaluated.

    A, b = SSpj_getFinal_A_b(M_C_K_Array, SolutionArray, f_Array, Δt, θ_Array)
"""
function SSpj_getFinal_A_b(M_C_K_Array::Array{SparseMatrixCSC{Float64,Int64},1}, SolutionArray::Array{Array{Float64, 1},1},
    f_Array::Array{Array{Float64,1},1}, Δt::Float64, θ_Array::Array{Float64,1})

    A::SparseMatrixCSC, SolMeanArray::Array{Array{Float64, 1},1}, f_mean::Array{Float64, 1} = get_SSpj_A_meanU_f(
    M_C_K_Array, SolutionArray, f_Array, Δt, θ_Array)
    lnM_C_K_Array::Int64 = length(M_C_K_Array)
    lnSolMeanArray::Int64 = length(SolutionArray)
    for i ∈ 1:lnSolMeanArray
        f_mean += M_C_K_Array[lnM_C_K_Array]*SolMeanArray[i]
        #println("f_mean at i: ", i, " =  ", f_mean[1])
        lnM_C_K_Array -=1
    end
    return -A, f_mean
end

"""Updates solutions: u, u̇, ü, etc; after obtaining the solution of α

    updateSolution!(SolutionArray::Array{Array{Float64, 1},1}, Δt::Float64, α::Array{Float64, 1})
"""
function updateSolution!(SolutionArray::Array{Array{Float64, 1},1}, Δt::Float64, α::Array{Float64, 1})

    for dot ∈ 1:length(SolutionArray)
        p::Int64 = 1
        for i ∈ dot+1:length(SolutionArray)#(dot-1)
            SolutionArray[dot] += Δt^p/factorial(p)*SolutionArray[i]
            p +=1
        end
        SolutionArray[dot] += Δt^p/factorial(p)*α
    end
    return nothing
end


"""Updates array of f vectors after obtaining the solution of f at a new time step.

    update_f!(f_Array::Array{Array{Float64,1},1}, f_new::Array{Float64,1})
"""
function update_f!(f_Array::Array{Array{Float64,1},1}, f_new::Array{Float64,1})
    popfirst!(f_Array)
    push!(f_new)
    return nothing
end
