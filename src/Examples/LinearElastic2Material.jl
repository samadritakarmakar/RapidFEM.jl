using RapidFEM, SparseArrays, WriteVTK

function LinearElastic2Material()
    mesh::Mesh = RapidFEM.readMesh("../test/Bar2.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib1::Tuple{Int64, Int64} = (3,4)
    volAttrib2::Tuple{Int64, Int64} = (3,5)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E1::Float64 = 200e3 #MPa
    ν1::Float64 = 0.3
    E2::Float64 = 1000e3 #MPa
    ν2::Float64 = 0.3
    tensorMap::Dict{Int64, Int64} = RapidFEM.getTensorMapping()
    C1::Array{Float64,2} = RapidFEM.createVoigtElasticTensor(E1, ν1)
    C2::Array{Float64,2} = RapidFEM.createVoigtElasticTensor(E2, ν2)
    K1::SparseMatrixCSC = RapidFEM.assembleMatrix((tensorMap, C1), volAttrib1, FeSpace, mesh, RapidFEM.local_∇v_C_∇u, problemDim, activeDimensions)
    K2::SparseMatrixCSC = RapidFEM.assembleMatrix((tensorMap, C2), volAttrib2, FeSpace, mesh, RapidFEM.local_∇v_C_∇u, problemDim, activeDimensions)
    K::SparseMatrixCSC = K1+K2
    source(x) = [0.0, 0.0, 0.0]
    f1::Vector = RapidFEM.assembleVector(source, volAttrib1, FeSpace, mesh, RapidFEM.localSource, problemDim, activeDimensions)
    f::Vector = f1 + RapidFEM.assembleVector(source, volAttrib2, FeSpace, mesh, RapidFEM.localSource, problemDim, activeDimensions)
    neumann(x) = [0.0, 0.0, -0.3333333333]
    f += RapidFEM.assembleVector(neumann, neumAttrib, FeSpace, mesh, RapidFEM.localNeumann, problemDim, activeDimensions)
    DirichletFunction(x) = zeros(problemDim)
    RapidFEM.applyDirichletBC!(K, f, DirichletFunction, dirchAttrib, mesh, problemDim)
    x::Vector = K\f
    vtkfile = RapidFEM.InitializeVTK(x, "LinearElastic2Material",mesh, [volAttrib1, volAttrib2], problemDim)
    vtkfile["Displacement"] = x
    RapidFEM.vtkSave(vtkfile)
    return nothing
end
