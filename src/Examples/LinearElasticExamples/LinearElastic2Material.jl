#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
using RapidFEM, SparseArrays, WriteVTK


function LinearElastic2Material()
    mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/Bar2.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib1::Tuple{Int64, Int64} = (3,4)
    volAttrib2::Tuple{Int64, Int64} = (3,5)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E1::Float64 = 200e3 #MPa
    ν1::Float64 = 0.3
    E2::Float64 = 300e3 #MPa
    ν2::Float64 = 0.3
    tensorMap::Dict{Int64, Int64} = RapidFEM.getTensorMapping()
    C1::Array{Float64,2} = RapidFEM.createVoigtElasticTensor(E1, ν1)
    C2::Array{Float64,2} = RapidFEM.createVoigtElasticTensor(E2, ν2)
    K::SparseMatrixCSC = RapidFEM.assembleMatrix((tensorMap, C1), volAttrib1, FeSpace, mesh, RapidFEM.local_∇v_C_∇u!, problemDim, activeDimensions)
    K += RapidFEM.assembleMatrix((tensorMap, C2), volAttrib2, FeSpace, mesh, RapidFEM.local_∇v_C_∇u!, problemDim, activeDimensions)
    source(x, varArgs...) = [0.0, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib1, FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    f += RapidFEM.assembleVector(source, volAttrib2, FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    neumann(x, varArgs...) = [0.0, 0.0, -0.3333333333]
    f += RapidFEM.assembleVector(neumann, neumAttrib, FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
    DirichletFunction(x, varArgs...) = zeros(problemDim)
    K = RapidFEM.applyDirichletBC!(f, K, DirichletFunction, dirchAttrib, mesh, problemDim)
    x::Vector = K\f
    σTemp::Array{Float64,1} = RapidFEM.InvDistInterpolation([RapidFEM.gaussianStress, RapidFEM.gaussianStress],
    x, [(tensorMap, C1), (tensorMap, C1)],  FeSpace, mesh,
    [volAttrib1, volAttrib2], problemDim, activeDimensions)
    σ::Array{Float64,1} = RapidFEM.voigtToTensor(σTemp, mesh)
    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("LinearElastic2Material", mesh, [volAttrib1, volAttrib2], problemDim)
    RapidFEM.vtkDataAdd!(vtkMeshData, (x, σ), ("Displacement", "Stress"))
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
end
