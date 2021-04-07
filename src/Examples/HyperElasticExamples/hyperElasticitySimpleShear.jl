using RapidFEM, LargeDefs, SparseArrays, LinearAlgebra, Tensors, NLsolve#, PyPlot
include("../hyperElasticLocalAssembly/hyperElastic.jl")

function hyperElasticity()
    #mesh::Mesh = RapidFEM.readMesh("../../test/OneElmntMsh/HexahedralOrder2.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cubeTest.msh")
    mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cube1elmnt.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/wedge.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,7)

    dirchAttribx::Tuple{Int64, Int64} = (2,1) #Lockx # Move this!!!!
    dirchAttribx_move::Tuple{Int64, Int64} = (2,2) #moveSurf # Move this!!!!
    dirchAttriby::Tuple{Int64, Int64} = (2,3) #Locky # Move this!!!!
    dirchAttriby_lock::Tuple{Int64, Int64} = (2,4) #Move_y #Lock this!!!
    dirchAttribz::Tuple{Int64, Int64} = (2,5) #locked_z #Lock this!!!
    dirchAttribz_lock::Tuple{Int64, Int64} = (2,6) #move_z #Lock this!!!
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 10.0 #MPa
    ν::Float64 = 0.3
    #λ = (ν*E)/((1+ν)*(1-2*ν))
    #μ = E/(2*(1+ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))

    println("λ = ", λ, " μ = ", μ)


    Dx = 0.5
    noOfSteps = 5

    b = [0.0, 0.0, 0.0]
    γ = zeros(noOfSteps)


    #hyperModel = LargeDefs.saintVenant
    hyperModel = LargeDefs.neoHookeanCompressible
    modelParams::Tuple  = (λ, μ)

    totalDoF::Int64 = mesh.noOfNodes*problemDim
    initSoln::Array{Float64, 1} = zeros(totalDoF)

    DirichletFunction1(x; varArgs...) = zeros(problemDim)

    assemble_f(initSoln) = begin
        start = time()
        source(x, varArgs...) = [b...]
         f::Array{Float64,1} = -RapidFEM.assembleVector((source, initSoln),
         volAttrib, FeSpace, mesh, localReferenceSource!, problemDim, activeDimensions)
        #println("Fx = ", Fx)

        hyperElasticParameters = (hyperModel, modelParams, initSoln)
         fσ::Array{Float64,1} = RapidFEM.assembleVector!(hyperElasticParameters,
         volAttrib, FeSpace, mesh, local_δE_S_Vector!, problemDim, activeDimensions)
        #println("\nnorm(fσ) = ", fσ,"\n")#, "\nnorm(initSoln) = ", initSoln)

         f+= fσ
         RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttribx, mesh,
         problemDim, [1,1,1])
         RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttribx_move, mesh,
         problemDim, [1,1,1])
         RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttriby, mesh,
         problemDim, [1,1,1])
         RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttriby_lock,
         mesh, problemDim, [1, 1, 1])
         RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttribz,
         mesh, problemDim, [0, 0, 1])
         RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttribz_lock,
         mesh, problemDim, [0, 0, 1])
        elasped = time() - start
        #println("Time for vector Assembly = ", elasped)
        return f
    end

    assemble_J(initSoln) = begin
        start = time()
        hyperElasticParameters = (hyperModel, modelParams, initSoln)
         J::SparseMatrixCSC = RapidFEM.assembleMatrix!(hyperElasticParameters,
         volAttrib, FeSpace, mesh, local_δE_Cᵀ_ΔE!, problemDim, activeDimensions)

        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttribx,
        mesh, problemDim, [1, 1, 1])
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttribx_move,
        mesh, problemDim, [1, 1, 1])
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttriby,
        mesh, problemDim, [1, 1, 1])
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttriby_lock,
        mesh, problemDim, [1, 1, 1])
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttribz,
        mesh, problemDim, [0, 0, 1])
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttribz_lock,
        mesh, problemDim, [0, 0, 1])
        elasped = time() - start
        #println("J = ", Matrix(J))
        #println("Time for Matrix Assembly = ", elasped)
        return J
    end

    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("cubeTestShear", mesh, [volAttrib], problemDim)
    #plot()
    for i ∈ 1:noOfSteps

        γ[i] = (i/noOfSteps)#^0.45

        println("i = ", i, " Dx = ", Dx)
        DirichletFunction_x(x; varArgs...) = [γ[i]*x[2]*Dx, 0.0, 0.0]
        DirichletFunction_Top(x; varArgs...) = [γ[i]*Dx, 0.0, 0.0]

        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction_x,
        dirchAttribx, mesh, problemDim, [1,1,1])
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction_x,
        dirchAttribx_move, mesh, problemDim, [1,1,1])
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction_Top,
        dirchAttriby, mesh, problemDim, [1,1,1])
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction1,
        dirchAttriby_lock, mesh, problemDim, [1,1,1])
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction1,
        dirchAttribz, mesh, problemDim, [0,0,1])
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction1,
        dirchAttribz_lock, mesh, problemDim, [0,0,1])

         #initSoln, convergenceData = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
        #    xtol = 1e-11, ftol = 1.e-5, relTol= 1e-8, iterations = 100, skipJacobian = 1 , printConvergence = true)

        solverResults = nlsolve(assemble_f, assemble_J, initSoln, show_trace = true,
        method = :trust_region, linesearch = BackTracking())
        initSoln .= solverResults.zero

        hyperElasticParameters = (hyperModel, modelParams)

        SecondPiolaStress = RapidFEM.InvDistInterpolation([gaussianSecondPiolaStress],
        initSoln, [hyperElasticParameters], FeSpace, mesh, [volAttrib],
            problemDim, activeDimensions)
        qsnDim = 1
        #for i ∈ 1:length(mesh.Elements[volAttrib])
        #    println("S[$qsnDim] = ", SecondPiolaStress[9*(i-1)+qsnDim])
        #end
        DeformationGrad = RapidFEM.InvDistInterpolation([gaussianDeformationGrad],
         initSoln, [hyperElasticParameters], FeSpace, mesh, [volAttrib],
             problemDim, activeDimensions)
        GreenStrain = RapidFEM.InvDistInterpolation([gaussianGreenStrain],
        initSoln, [hyperElasticParameters], FeSpace, mesh, [volAttrib],
                problemDim, activeDimensions)
        RapidFEM.vtkDataAdd!(vtkMeshData, (initSoln, SecondPiolaStress,
        GreenStrain, DeformationGrad), ("Displacement", "Second Piola Stress",
        "Green Strain", "Deformation Gradient"), float(i), i)
        #plot(log10.(convergenceData.relNorm), linestyle= :dashdot)
    end
    #xlabel("Iterations")
    #ylabel("Log base 10 of Relative Convergence")
    RapidFEM.vtkSave(vtkMeshData)
end
