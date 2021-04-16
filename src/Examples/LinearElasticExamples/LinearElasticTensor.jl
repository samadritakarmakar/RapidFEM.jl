#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

using RapidFEM, SparseArrays, WriteVTK, LinearAlgebra, IterativeSolvers, Tensors

include("../LinearElasticLocalAssembly/linearElasticity.jl")

function LinearElastic()
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/BarFine.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cubeTestO1.msh")
    mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cubeTestPntO1.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cubeTestR1.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,5)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 1 #MPa
    ν::Float64 = 0.3
    C = createElasticTensor(E, ν)
    K::SparseMatrixCSC = RapidFEM.assembleMatrix!((C, ), volAttrib, FeSpace, mesh, local_∇v_C_∇u_Tensor!, problemDim, activeDimensions)
    #println(Matrix(K))
    stiffStringArray = string.(Matrix(K))
    writeData = ""
    for i ∈ 1:size(stiffStringArray,1)
      for j ∈ 1:size(stiffStringArray,2)
        writeData = writeData*stiffStringArray[i,j]*"  "
      end
      writeData = writeData*"\n"
    end
    file = open("stiffness", "w")
    write(file, writeData)
    close(file)
    source(x, varArgs...) = [0.0, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector!(source, volAttrib, FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    neumann(x; varArgs...) = [0.0, 0.0, 0.0]
    f += RapidFEM.assembleVector!(neumann, neumAttrib, FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
    DirichletFunction(x; varArgs...) = zeros(problemDim)
    RapidFEM.applyDirichletBC!(f, K, DirichletFunction, dirchAttrib, mesh, problemDim)
    DirichletFunction2(x; varArgs...) = begin
      if (abs(x[1]-1.0)<1e-14 && abs(x[2]-1.0)<1e-14 && abs(x[3]-1.0)<1e-14)
        return [0.1, 0.1, 0.1]
      end
      return zeros(3)
    end
    RapidFEM.applyDirichletBC!(f, K, DirichletFunction2, (0,6), mesh, problemDim)
    x::Vector = K\f
    println(x)
    #x = cg(K,f)
    σ::Array{Float64,1} = RapidFEM.InvDistInterpolation([gaussianStress],
    x, [(C,)],  FeSpace, mesh,  [volAttrib],
    problemDim, activeDimensions, kwargs=0)
    ϵ::Array{Float64,1} = RapidFEM.InvDistInterpolation([gaussianStrain],
    x, [(C,)],  FeSpace, mesh,  [volAttrib],
    problemDim, activeDimensions, kwargs=0)

    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("LinearElasticR0O1", mesh, [volAttrib], problemDim)
    RapidFEM.vtkDataAdd!(vtkMeshData, (x, σ, ϵ), ("Displacement", "Stress","Strain"))
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
end
