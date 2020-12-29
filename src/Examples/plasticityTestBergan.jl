#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
using RapidFEM, SmallStrainPlastic, SparseArrays, WriteVTK, NLsolve, LinearAlgebra, Plots
include("plasticLocalAssembly/smallStrainPlasticity.jl")

function plasticity()
    #mesh::Mesh = RapidFEM.readMesh("../test/MeshFiles/Bar.msh")
    mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/HexahedralL3O2.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,5)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttribx::Tuple{Int64, Int64} = (2,1) #Lockx
    dirchAttriby::Tuple{Int64, Int64} = (2,3) #Locky
    dirchAttribz::Tuple{Int64, Int64} = (2,4) #Lockz
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    σ_y::Float64 = 200.0
    Fx₀::Float64 = 400.0
    maxLoadLimit = 210.0
    minLoadLimit = 0.0

    forceArray = zeros(0)

    residualArray::Array{Float64, 1} = zeros(0)
    tensorMap::Dict{Int64, Int64} = RapidFEM.getTensorMapping()
    C::Array{Float64, 2} = SmallStrainPlastic.getMandelElasticTensor(E, ν)

    #Intializing SmallStrainPlastic Library
    model::PlasticModel = SmallStrainPlastic.j2Model
    stateDict = SmallStrainPlastic.createStateDict()
    stateDictBuffer = SmallStrainPlastic.createStateDict()
    params_J2 = SmallStrainPlastic.initParams_j2(σ_y, 20e3)

    #λ = 1.0
    Δλ::Float64 = 1.0

    totalDoF::Int64 = mesh.noOfNodes*problemDim
    initSoln::Array{Float64, 1} = zeros(totalDoF)
    finalSoln::Array{Float64, 1} = zeros(totalDoF)
    f₀::Array{Float64, 1} = zeros(totalDoF)
    fIter = 0
    mainIter = 1
    DirichletFunction(x; varArgs...) = zeros(problemDim)

    applyDirchletBC_on_f!(f) = begin
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttribx, mesh, problemDim, [1, 0, 0])
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttriby, mesh, problemDim, [0, 1, 0])
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttribz, mesh, problemDim, [0, 0, 1])
        return nothing
    end

    applyDirchletBC_onSoln!(initSoln) = begin
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttribx,
        mesh, problemDim, [1,0,0])
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttriby,
        mesh, problemDim, [0,1,0])
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttribz,
        mesh, problemDim, [0,0,1])
    end

    assemble_f₀() = begin
        source(x, varArgs...) = [0.0, 0.0, 0.0]
        f::Array{Float64,1} = -RapidFEM.assembleVector(source, volAttrib, FeSpace,
        mesh, RapidFEM.localSource!, problemDim, activeDimensions)
        println("Fx = ", Fx₀)
        neumann(x; varArgs...) = [Fx₀, 0.0, 0.0]
        f -= RapidFEM.assembleVector(neumann, neumAttrib, FeSpace,
        mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
        applyDirchletBC_on_f!(f)
        return f
    end

    assemble_f(initSoln) = begin

        tensorMap_N_PlasticData = (tensorMap, C, model, params_J2, stateDict, stateDictBuffer, initSoln)
        fσ::Array{Float64,1} = RapidFEM.assembleVector!(tensorMap_N_PlasticData, volAttrib, FeSpace,
        mesh, local_∇v_σ_Vector!, problemDim, activeDimensions)
        applyDirchletBC_on_f!(fσ)

        if (fIter > 0)
            Δλ = RapidFEM.berganIncrement(f₀, fσ)
        end
        println("Δλ =", Δλ, " norm(fσ) = ", norm(fσ), " norm(Δλ*f₀) = ", norm(Δλ*f₀))
        f = fσ + Δλ*f₀

        fIter += 1
        return f
    end

    assemble_J(initSoln) = begin
        tensorMap_N_PlasticData = (tensorMap, C, model, params_J2, stateDict, stateDictBuffer, initSoln)
        J::SparseMatrixCSC = RapidFEM.assembleMatrix!(tensorMap_N_PlasticData, volAttrib, FeSpace,
        mesh, local_∇v_Cᵀ_∇u!, problemDim, activeDimensions)
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttribx, mesh, problemDim, [1,0,0])
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttriby, mesh, problemDim, [0,1,0])
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttribz, mesh, problemDim, [0,0,1])
        #println(J)
        return J
    end

    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("Plasticity", mesh, [volAttrib], problemDim)

    writeToParaview(initSoln, mainIter) = begin
        tensorMap_N_PlasticData = (tensorMap, C, model, params_J2, stateDict, stateDictBuffer, initSoln)
        ϵᵖTemp::Array{Float64,1} = RapidFEM.InvDistInterpolation([gaussian_ϵᵖ],
        initSoln, [tensorMap_N_PlasticData],  FeSpace, mesh,  [volAttrib],
        problemDim, activeDimensions)
        ϵᵖ::Array{Float64,1} = RapidFEM.voigtToTensor(ϵᵖTemp, mesh)

        ϵ::Array{Float64,1} = RapidFEM.InvDistInterpolation([gaussian_ϵ],
        initSoln, [tensorMap_N_PlasticData],  FeSpace, mesh,  [volAttrib],
        problemDim, activeDimensions)
        #ϵ::Array{Float64,1} = RapidFEM.voigtToTensor(ϵTemp, mesh)
        #println("finalSoln = ", finalSoln)
        σ::Array{Float64,1} = RapidFEM.InvDistInterpolation([gaussian_σ],
        initSoln, [tensorMap_N_PlasticData],  FeSpace, mesh,  [volAttrib],
        problemDim, activeDimensions)
        #σ::Array{Float64,1} = RapidFEM.voigtToTensor(σTemp, mesh)
        println("σ = ", σ[1])
        RapidFEM.vtkDataAdd!(vtkMeshData, (initSoln,ϵᵖ, ϵ, σ),
        ("Displacement", "PlasticStrain", "Strain", "Stress"), float(mainIter),mainIter)
    end


    for cycle ∈ 1:2
        f₀ = assemble_f₀()
        while(abs(Δλ-1.0)>1e-9 || mainIter < 3)
            println("\nmainIter = ", mainIter, " abs(Δλ-1.0) = ", abs(Δλ-1.0))
            Δλ =1.0
            applyDirchletBC_onSoln!(initSoln)
            #solverResults = NLsolve.nlsolve(assemble_f, assemble_J, initSoln;
            #xtol = 1e-7, ftol = 1e-7, iterations = 1000)
            initSoln = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
                xtol = 1e-12, ftol = 1e-8, iterations = 2000, skipJacobian =1 , printConvergence = true)
            fIter = 0
            mainIter +=1
            #initSoln .= finalSoln
            SmallStrainPlastic.updateStateDict4rmBuffer!(stateDict, stateDictBuffer)
            writeToParaview(initSoln, mainIter)
        end
        Fx₀ = 1e-8
        Δλ = 0.5
    end

    RapidFEM.vtkSave(vtkMeshData)
    #return nothing
    #plot(forceArray, label = ["Force"])

end
