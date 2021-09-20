#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

using RapidFEM, SparseArrays, WriteVTK, LinearAlgebra, IterativeSolvers
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

function LinearElastic()
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/Bar.msh")
    mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/BarNew1_R2_O2.msh")
    
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cubeTestO1.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cubeTestPntO1.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,4)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    tensorMap::Dict{Int64, Int64} = RapidFEM.getTensorMapping()
    C::Array{Float64,2} = RapidFEM.createVoigtElasticTensor(E, ν)
    K::SparseMatrixCSC = RapidFEM.assembleMatrix((tensorMap, C), volAttrib, FeSpace, mesh, RapidFEM.local_∇v_C_∇u!, problemDim, activeDimensions)
    source(x, varArgs...) = [0.0, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib, FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    neumann(x; varArgs...) = [0.0, 0.0, 0.333333333]
    f += RapidFEM.assembleVector(neumann, neumAttrib, FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
    DirichletFunction(x; varArgs...) = zeros(problemDim)
    K = RapidFEM.applyDirichletBC!(f, K, DirichletFunction, dirchAttrib, mesh, problemDim)
    println(size(K))
    x::Vector = K\f
    #println(x)
    #x = cg(K,f)
    σTemp::Array{Float64,1} = RapidFEM.InvDistInterpolation([RapidFEM.gaussianStress],
    x, [(tensorMap, C)],  FeSpace, mesh,  [volAttrib],
    problemDim, activeDimensions, kwargs=0)
    σ::Array{Float64,1} = RapidFEM.voigtToTensor(σTemp, mesh)
    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("LinearElastic",mesh, [volAttrib], problemDim)
    RapidFEM.vtkDataAdd!(vtkMeshData, (x, σ), ("Displacement", "Stress"))
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
end
