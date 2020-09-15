using RapidFEM, SparseArrays, WriteVTK, LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

function LinearElastic()
    mesh::Mesh = RapidFEM.readMesh("../test/BarHex.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    tensorMap::Dict{Int64, Int64} = RapidFEM.getTensorMapping()
    C::Array{Float64,2} = RapidFEM.createVoigtElasticTensor(E, ν)
    K::SparseMatrixCSC = RapidFEM.assembleMatrix((tensorMap, C), volAttrib, FeSpace, mesh, RapidFEM.local_∇v_C_∇u, problemDim, activeDimensions)
    source(x) = [0.0, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib, FeSpace, mesh, RapidFEM.localSource, problemDim, activeDimensions)
    neumann(x) = [0.0, 0.0, -0.3333333333]
    f += RapidFEM.assembleVector(neumann, neumAttrib, FeSpace, mesh, RapidFEM.localNeumann, problemDim, activeDimensions)
    DirichletFunction(x) = zeros(problemDim)
    RapidFEM.applyDirichletBC!(K, f, DirichletFunction, dirchAttrib, mesh, problemDim)
    x::Vector = K\f
    vtkfile = RapidFEM.InitializeVTK(x, "LinearElastic",mesh, [volAttrib], problemDim)
    vtkfile["Displacement"] = x
    σTemp::Array{Float64,1} = RapidFEM.InvDistInterpolation(RapidFEM.gaussianStress, x, (tensorMap, C),  FeSpace, mesh,  volAttrib, problemDim, activeDimensions)
    σ::Array{Float64,1} = RapidFEM.voigtToTensor(σTemp, mesh)
    vtkfile["Stress"] = σ
    RapidFEM.vtkSave(vtkfile)
    return nothing
end
