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

"""Scales the mode shapes after, Normalizing and scaling by the scalingConstant, each of the mode shapes. 
Normalization with scaling is done by, uᴿ = scalingConstant*û/u_L2Norm, where u_L2Norm = √(∫ û⋅û dΩ)"""
function scaleModeShapes!(u_modeShapes::AbstractMatrix{Float64}, volAttrib::Tuple, FeSpace::Dict, mesh::Mesh, 
  problemDim::Int64, activeDims::Vector{Int64}, scalingConstant::Number = 1.0)
  
  dimRange = RapidFEM.createDimRange()
  usedDims = dimRange[activeDims]
  modeNo = 1
  for modeShape ∈ eachcol(u_modeShapes)
    elementNo = 1
    L2Norm = 0.0
    for element ∈ mesh.Elements[volAttrib]
      coordArray = (getCoordArray(mesh, element))[usedDims, :]
      shapeFunction = feSpace!(FeSpace, element, mesh)
      modeAtElNodes = getSolAtElement(modeShape, element, problemDim, activeDims)
      #ipNo = 1
      noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
      modeAtIp = zeros(problemDim)
      for ipNo ∈ 1:noOfIpPoints
          ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
          dΩ = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
          ϕ = get_ϕ(shapeFunction, ipNo)
          get_u!(modeAtIp, modeAtElNodes, ϕ, problemDim)
          L2Norm += dot(modeAtIp, modeAtIp)*dΩ
          #ipNo += 1
      end
      elementNo += 1
    end
    L2Norm = sqrt(L2Norm)
    #println("L2Norm of mode no $modeNo = $L2Norm")
    u_modeShapes[:, modeNo] .= scalingConstant/L2Norm .*modeShape
    modeNo += 1
  end
  return u_modeShapes
end
