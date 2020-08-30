using RapidFEM, SparseArrays, WriteVTK, FEMSparse

function poissonEquation()
    mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/HexahedralOrder1.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 1
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,2)
    dirchAttrib::Tuple{Int64, Int64} = (2,1)
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    parameterFunction(x) = [1.0]#, 1.0, 1.0]
    K::SparseMatrixCSC = RapidFEM.assembleMatrix(parameterFunction, volAttrib,
    FeSpace, mesh, RapidFEM.local_lagrange_K, problemDim, activeDimensions)
    source(x) = [0.0]#, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib,
    FeSpace, mesh, RapidFEM.localSource, problemDim, activeDimensions)
    neumann(x) = [0.1]#, 0.0, 0.0]
    f += RapidFEM.assembleVector(neumann, neumAttrib,
    FeSpace, mesh, RapidFEM.localNeumann, problemDim, activeDimensions)
    DirichletFunction(x) = zeros(problemDim)
    RapidFEM.applyDirichletBC!(K, f, DirichletFunction, dirchAttrib,
    mesh, problemDim)
    x::Vector = K\f
    vtkfile = RapidFEM.InitializeVTK(x, "poisson",mesh, [volAttrib],problemDim)
    vtkfile["Displacement"] = x
    RapidFEM.vtkSave(vtkfile)
    return FeSpace, K_local
end
