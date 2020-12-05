#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
using RapidFEM, SmallStrainPlastic, SparseArrays, WriteVTK, NLsolve, LinearAlgebra
include("plasticLocalAssembly/smallStrainPlasticity.jl")


function plasticity()
    #mesh::Mesh = RapidFEM.readMesh("../test/MeshFiles/Bar.msh")
    mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/TetrahedralOrder2.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    σ_y::Float64 = 200.0
    Fx::Float64 = 40.0
    residualArray::Array{Float64, 1} = zeros(0)
    tensorMap::Dict{Int64, Int64} = RapidFEM.getTensorMapping()
    C::Array{Float64, 2} = SmallStrainPlastic.createVoigtElasticTensor(E, ν)

    #Intializing SmallStrainPlastic Library
    model::PlasticModel = SmallStrainPlastic.j2Model
    stateDict = SmallStrainPlastic.createStateDict()
    stateDictBuffer = SmallStrainPlastic.createStateDict()
    params_J2 = SmallStrainPlastic.initParams_j2(σ_y, 20e3)

    totalDoF::Int64 = mesh.noOfNodes*problemDim
    #f::Array{Float64,1} = zeros(totalDoF)
    initSoln::Array{Float64, 1} = zeros(totalDoF)
    finalSoln::Array{Float64, 1} = zeros(totalDoF)
    #J::SparseMatrixCSC = spzeros(totalDoF, totalDoF)
    DirichletFunction(x; varArgs...) = zeros(problemDim)
#=
    assemble_fJ(initSoln) = begin
        source(x, varArgs...) = [0.0, 0.0, 0.0]
        f::Array{Float64,1} = -RapidFEM.assembleVector(source, volAttrib, FeSpace,
        mesh, RapidFEM.localSource!, problemDim, activeDimensions)
        #println("Fx = ", Fx)
        neumann(x; varArgs...) = [Fx, 0.0, 0.0]
        f -= RapidFEM.assembleVector(neumann, neumAttrib, FeSpace,
        mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
        #println("external = ", f)

        tensorMap_N_PlasticData = (tensorMap, plasticVars, RapidFEM.j2Model, params_J2, initSoln)

        fσ::Array{Float64,1}, J::SparseMatrixCSC = RapidFEM.assembleVectorMatrix!(tensorMap_N_PlasticData, volAttrib, FeSpace,
        mesh, RapidFEM.local_∇v_σ_Vector!, RapidFEM.local_∇v_Cᵀ_∇u!, problemDim, activeDimensions)

        f+= fσ

        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttrib, mesh, problemDim)
        J = applyNLDirichletBC_on_J!(J, dirchAttrib, mesh, problemDim)

        println("norm(f) = ", norm(f))#, "\nnorm(fσ) = ", fσ, "\nnorm(initSoln) = ", initSoln)
        #println("\nnorm(initSoln) = ", initSoln)
        finalSoln .= initSoln
        #println("norm(J) = ", Matrix(J))
        return f, J
    end=#

    assemble_f(initSoln) = begin
        source(x, varArgs...) = [0.0, 0.0, 0.0]
        f::Array{Float64,1} = -RapidFEM.assembleVector(source, volAttrib, FeSpace,
        mesh, RapidFEM.localSource!, problemDim, activeDimensions)
        #println("Fx = ", Fx)
        neumann(x; varArgs...) = [Fx, 0.0, 0.0]
        f -= RapidFEM.assembleVector(neumann, neumAttrib, FeSpace,
        mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
        println("external = ", f)

        tensorMap_N_PlasticData = (tensorMap, C, model, params_J2, stateDict, stateDictBuffer, initSoln)
        fσ::Array{Float64,1} = RapidFEM.assembleVector!(tensorMap_N_PlasticData, volAttrib, FeSpace,
        mesh, local_∇v_σ_Vector!, problemDim, activeDimensions)

        f+= fσ
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttrib, mesh, problemDim)
        #RapidFEM.applyNLDirichletBC_on_f!(fσ, dirchAttrib, mesh, problemDim)
        println("\nnorm(fσ) = ", fσ,"\n")#, "\nnorm(initSoln) = ", initSoln)
        return f
    end

    assemble_J(initSoln) = begin
        tensorMap_N_PlasticData = (tensorMap, C, model, params_J2, stateDict, stateDictBuffer, initSoln)
        J::SparseMatrixCSC = RapidFEM.assembleMatrix!(tensorMap_N_PlasticData, volAttrib, FeSpace,
        mesh, local_∇v_Cᵀ_∇u!, problemDim, activeDimensions)
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttrib, mesh, problemDim)
        return J
    end

    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("Plasticity", mesh, [volAttrib], problemDim)
    for i ∈ 1:20
        println("i = ", i)
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttrib,
        mesh, problemDim)
        #solverResults = NLsolve.nlsolve(assemble_f, assemble_J, initSoln;
        #xtol = 1e-7, ftol = 1e-7, iterations = 30)
        finalSoln = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
            xtol = 1e-7, ftol = 1e-7, iterations = 50, printConvergence = true)
        initSoln .= finalSoln

        SmallStrainPlastic.updateStateDict4rmBuffer!(stateDict, stateDictBuffer)
        #println("finalSoln = ", finalSoln)
        RapidFEM.vtkDataAdd!(vtkMeshData, (finalSoln,),
        ("Displacement", ), float(i), i)
        if i<=15
            Fx = 40*(i+1.0)
        else
            Fx = 0.0
        end
    end
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
    #plot(residualArray)

end
