using RapidFEM, DifferentialEquations

function assembleMassMatrix!(FeSpace, mesh::Mesh, domainAttrib::Tuple{Int64, Int64},
    densityFunction::Function, problemDim::Int64 = 3, activeDimensions::Array{Int64,1} = [1, 1, 1])

    M::SparseMatrixCSC = RapidFEM.assembleMatrix(densityFunction, volAttrib,
    FeSpace, mesh, RapidFEM.local_v_ρ_u!, problemDim, activeDimensions)

    return M
end

function assembleStiffMatrix!(FeSpace, mesh::Mesh, domainAttrib::Tuple{Int64, Int64},
    conductivityFunction::Function, problemDim::Int64 = 3, activeDimensions::Array{Int64,1} = [1, 1, 1])

    K::SparseMatrixCSC = RapidFEM.assembleMatrix(conductivityFunction, volAttrib,
    FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)

    return K
end

function assembleSurfTrnsfrMatrix!(FeSpace, mesh::Mesh, domainAttrib::Tuple{Int64, Int64},
    htTrnsfrCoeffFunc::Function, problemDim::Int64 = 3, activeDimensions::Array{Int64,1} = [1, 1, 1])

    K2::SparseMatrixCSC = RapidFEM.assembleMatrix(htTrnsfrCoeffFunc, volAttrib,
    FeSpace, mesh, RapidFEM.localBoundary_v_ρ_u!, problemDim, activeDimensions)

    return K2
end

function assembleSourceVector!(FeSpace, mesh::Mesh, domainAttrib::Tuple{Int64, Int64},
    sourceFunc::Function, problemDim::Int64 = 3, activeDimensions::Array{Int64,1} = [1, 1, 1])

    f::Vector = RapidFEM.assembleVector(sourceFunc, volAttrib,
    FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)

    return f
end

function assembleNeumannVector!(FeSpace, mesh::Mesh, domainAttrib::Tuple{Int64, Int64},
    neumannFunc::Function, problemDim::Int64 = 3, activeDimensions::Array{Int64,1} = [1, 1, 1])

    f::Vector = RapidFEM.assembleVector(neumannFunc, volAttrib,
    FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)

    return f
end
