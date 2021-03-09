#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
using RapidFEM, SparseArrays, WriteVTK, FEMSparse

function poissonEquation()
    #mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/TetrahedralOrder2.msh")
    mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/HexahedralOrder1.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 1
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,2)
    dirchAttrib::Tuple{Int64, Int64} = (2,1)
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    parameterFunction(x, varArgs...) = [1.0]#, 1.0, 1.0]
    K::SparseMatrixCSC = RapidFEM.assembleMatrix(parameterFunction, volAttrib,
    FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)
    source(x, varArgs...) = [0.0]#, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib,
    FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    neumann(x, varArgs...) = [0.1]#, 0.0, 0.0]
    f += RapidFEM.assembleVector(neumann, neumAttrib,
    FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
    DirichletFunction(x, varArgs...) = zeros(problemDim)
    K = RapidFEM.applyDirichletBC!(f, K, DirichletFunction, dirchAttrib,
    mesh, problemDim)
    x::Vector = K\f
    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("poisson", mesh, [volAttrib],problemDim)
    RapidFEM.vtkDataAdd!(vtkMeshData, (x,), ("Field",))
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
end
