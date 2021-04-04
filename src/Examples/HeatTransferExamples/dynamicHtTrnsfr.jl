#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
using RapidFEM, SparseArrays

function dynHeatTransfer()
    mesh::Mesh = RapidFEM.readMesh("../test/MeshFiles/Rod.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 1
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,2)
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    #Mass Matrix
    ρ_C_Function(x, varArgs...) = [1.0]
    C::SparseMatrixCSC = RapidFEM.assembleMatrix(ρ_C_Function, volAttrib,
        FeSpace, mesh, RapidFEM.local_v_ρ_u!, problemDim, activeDimensions)
    #Stifffness Matrix
    conductivityFunction(x, varArgs...) = [1.0]
    K::SparseMatrixCSC = RapidFEM.assembleMatrix(conductivityFunction, volAttrib,
        FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)
    #Surface heat transfer
    htTrnsfrCoeffFunc(x, varArgs...) = [0.0]
    K += RapidFEM.assembleMatrix(htTrnsfrCoeffFunc, neumAttrib,
        FeSpace, mesh, RapidFEM.localBoundary_v_ρ_u!, problemDim, activeDimensions)
    #source vector
    source(x, varArgs...) = [0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib,
        FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    #neumann vector
    neumann(x, varArgs...) = [0.1]
    f -= RapidFEM.assembleVector(neumann, neumAttrib,
        FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)

    #Initial BC condition
    InitialBcFunc(x, varArgs...) = [0.0]
    u = zeros(length(f))
    RapidFEM.applyInitialBC!(u, InitialBcFunc,
    volAttrib, mesh, problemDim)
    #Time Integration parameters
    SolutionArray = [u]
    f_Array = [f, f]
    Δt = .0001
    θ_Array = [0.878]
    #Intiailization of Paraview output
    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("dynamicHeatTransfer", mesh, [volAttrib], problemDim)
    RapidFEM.vtkDataAdd!(vtkMeshData, (SolutionArray[1],), ("Temperature",), 0.0*Δt, 1)
    #Time Intengration
    DirichletFunction(x, varArgs...) = [100.0]
    for i ∈ 1:200
        println("At step: ", i, " time: ", i*Δt)
        A_global, f_mean_global = SSpj_getFinal_A_b([C, K], SolutionArray, f_Array, Δt, θ_Array)
        A_global = RapidFEM.applyDynamicDirichletBC!(SolutionArray, f_mean_global,
        A_global, DirichletFunction, dirchAttrib, mesh, problemDim)
        α_global = A_global\f_mean_global
        updateSolution!(SolutionArray, Δt, α_global)
        RapidFEM.vtkDataAdd!(vtkMeshData, (SolutionArray[1],), ("Temperature",), i*Δt, i+1)
    end
    RapidFEM.vtkSave(vtkMeshData)
end
