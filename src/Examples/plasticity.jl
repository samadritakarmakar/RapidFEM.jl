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
    mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/HexahedralOrder2.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    σ_y::Float64 = 200.0
    Fx::Float64 = 4.0
    maxLoadLimit = 250.0
    minLoadLimit = 0.0
    stepsMaxLoad = 40.0
    stepsMinLoad = 80.0
    stepMatrix = [0.0 0.0 1.0
                stepsMaxLoad stepsMaxLoad^2 1.0
                stepsMinLoad stepsMinLoad^2 1.0]
    stepCoeff = stepMatrix\[0.0; maxLoadLimit; minLoadLimit]

    steps = stepsMinLoad

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

    assemble_f(initSoln) = begin
        source(x, varArgs...) = [0.0, 0.0, 0.0]
        f::Array{Float64,1} = -RapidFEM.assembleVector(source, volAttrib, FeSpace,
        mesh, RapidFEM.localSource!, problemDim, activeDimensions)
        println("Fx = ", Fx)
        neumann(x; varArgs...) = [Fx, 0.0, 0.0]
        f -= RapidFEM.assembleVector(neumann, neumAttrib, FeSpace,
        mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
        #println("external = ", f)

        tensorMap_N_PlasticData = (tensorMap, C, model, params_J2, stateDict, stateDictBuffer, initSoln)
        fσ::Array{Float64,1} = RapidFEM.assembleVector!(tensorMap_N_PlasticData, volAttrib, FeSpace,
        mesh, local_∇v_σ_Vector!, problemDim, activeDimensions)

        f+= fσ
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttrib, mesh, problemDim)
        #RapidFEM.applyNLDirichletBC_on_f!(fσ, dirchAttrib, mesh, problemDim)
        #println("\nnorm(fσ) = ", fσ,"\n")#, "\nnorm(initSoln) = ", initSoln)
        #finalSoln .= initSoln
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

    for i ∈ 1:Int64(steps)
        FxTemp = [float(i) float(i)^2 1.0]*stepCoeff
        Fx = FxTemp[1]

        println("i = ", i)
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttrib,
        mesh, problemDim)
        #solverResults = NLsolve.nlsolve(assemble_f, assemble_J, initSoln;
        #xtol = 1e-7, ftol = 1e-7, iterations = 1000)
        initSoln = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
            xtol = 1e-12, ftol = 1e-8, iterations = 2000, skipJacobian =1 , printConvergence = true)
        #initSoln .= finalSoln

        SmallStrainPlastic.updateStateDict4rmBuffer!(stateDict, stateDictBuffer)

        tensorMap_N_PlasticData = (tensorMap, C, model, params_J2, stateDict, stateDictBuffer, initSoln)
        ϵᵖTemp::Array{Float64,1} = RapidFEM.InvDistInterpolation([gaussian_ϵᵖ],
        initSoln, [tensorMap_N_PlasticData],  FeSpace, mesh,  [volAttrib],
        problemDim, activeDimensions)
        ϵᵖ::Array{Float64,1} = RapidFEM.voigtToTensor(ϵᵖTemp, mesh)
        #println("finalSoln = ", finalSoln)
        RapidFEM.vtkDataAdd!(vtkMeshData, (initSoln,ϵᵖ),
        ("Displacement", "PlasticStrain"), float(i), i)

    end
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
    #plot(residualArray)

end
