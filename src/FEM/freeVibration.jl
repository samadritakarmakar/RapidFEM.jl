 """Applies locked boundary condition to 'attribute'. Designed to solve free vibration problem 
 in an inverse sense, meaning (1/ω^2) ū = K⁻¹M ū"""
 function applyFreeVibrationBC!(M::SparseMatrixCSC, K::SparseMatrixCSC,
    attribute::Tuple{Int64, Int64}, mesh::Mesh,
    problemDim::Int64, appliedDof::Array{Int64, 1} = ones(Int, problemDim))

    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim, appliedDof)
    P::SparseMatrixCSC, Pzeros::SparseMatrixCSC = getPermutionMatrix(vNodes, mesh, problemDim)
    M .= P*M*P
    K .= P*K*P
    K .+= Pzeros
    return M, K
 end

"""Given a certain vector of angular momentums, ω and a matrix of amplitudes, 
it gives the displacement at a certain time t, for a certain mode number"""
 function getFreeVibrationDisplacement(ω::Vector{Float64}, amplitudes::Matrix{Float64}, modeNumber::Int64, t::Float64)
   return amplitudes[:, modeNumber]*cos(ω[modeNumber]*t)
 end

"""After applying free vibration boundary conditions, this solves the problem using Arpack.eigs which can be used with
sparse matrices. Returns angular frequency, ω and amplitudes, ū"""
 function solveFreeVibration(M::AbstractMatrix{Float64}, K::AbstractMatrix{Float64}, noOfModes = 15)

  λ, u_amp = Arpack.eigs(Symmetric(M), Symmetric(K), maxiter = 1000, which = :LR, nev = noOfModes)
  return (1.0./sqrt.(λ)), u_amp
end

"""After applying free vibration boundary conditions, this solves the problem using the default julia Lapack backend.
Should not be used with sparse matrices. Returns angular frequency, ω and amplitudes, ū"""
function solveFreeVibrationDense(M::AbstractMatrix{Float64}, K::AbstractMatrix{Float64}, noOfModes = 15)

  λ, u_amp = eigen(Matrix(M), Matrix(K))
  return (1.0./sqrt.(λ))[end:-1:(length(λ)-noOfModes+1)], u_amp[:, end:-1:(length(λ)-noOfModes+1)]
end

function scaleModeShapes!(u_modeShapes::AbstractMatrix{Float64}, scalingConstant::Number = 1.0)
  noOfModes = size(u_modeShapes, 2)
  
  for modeNo ∈ 1:noOfModes
    scalingFactor = scalingConstant/norm(u_modeShapes[:,modeNo])
    u_modeShapes[:, modeNo] .= scalingFactor .* u_modeShapes[:, modeNo]
  end
  return u_modeShapes
end
