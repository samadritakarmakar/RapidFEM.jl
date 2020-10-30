using RapidFEM, DifferentialEquations

function dynHeatTransfer()
    mesh::Mesh = RapidFEM.readMesh("../test/MeshFiles/Bar.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 1
    volAttrib::Tuple{Int64, Int64} = (3,4)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]

    K::SparseMatrixCSC = RapidFEM.assembleMatrix(densityFunction, domainAttrib,
    FeSpace, mesh, RapidFEM.local_v_ρ_u!, problemDim, activeDimensions)

    K += RapidFEM.assembleMatrix(conductivityFunction, domainAttrib,
    FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)

    K += RapidFEM.assembleMatrix(htTrnsfrCoeffFunc, SurfHtLossAttrib,
    FeSpace, mesh, RapidFEM.localBoundary_v_ρ_u!, problemDim, activeDimensions)

    f::Vector = RapidFEM.assembleVector(sourceFunc, domainAttrib,
    FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)

    f += RapidFEM.assembleVector(neumannFunc, SurfHtLossAttrib,
    FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
end
