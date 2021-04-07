using RapidFEM, LargeDefs, SparseArrays, LinearAlgebra, Tensors, NLsolve, LineSearches#, PyPlot
include("../hyperElasticLocalAssembly/hyperElastic.jl")

function hyperElasticity()
    #mesh::Mesh = RapidFEM.readMesh("../../test/OneElmntMsh/HexahedralOrder2.msh")
    mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cube2.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cube1elmnt.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/wedge.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,4)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib1::Tuple{Int64, Int64} = (2,1) #Lock
    #dirchAttrib2::Tuple{Int64, Int64} = (2,2) #move
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 10.0 #MPa
    ν::Float64 = 0.3
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    #λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    #μ = E / (2(1 + ν))

    Fx::Float64 = 0.0
    Fy::Float64 = 0.0
    Fz::Float64 = 0.0

    Fx_max::Float64 = 0.0
    Fy_max::Float64 = 0.5
    Fz_max::Float64 = 0.5

    disp_max(X) = begin
        x, y, z = X
        L = 1.0
        θ = π/3
        0.5*[0.0,
        L/2 - y + (y-L/2)*cos(θ) - (z-L/2)*sin(θ),
        L/2 - z + (y-L/2)*sin(θ) + (z-L/2)*cos(θ)]
    end
    noOfSteps = 1

    b_max = [0.0, 0.0, 0.0]

    b = zeros(problemDim)
    γ = zeros(noOfSteps+1)


    #hyperModel = LargeDefs.saintVenant
    hyperModel = LargeDefs.neoHookeanCompressible
    modelParams::Tuple  = (λ, μ)

    totalDoF::Int64 = mesh.noOfNodes*problemDim
    initSoln::Array{Float64, 1} = zeros(totalDoF)

    DirichletFunction1(x; varArgs...) = zeros(problemDim)

    assemble_f(initSoln) = begin
        start = time()
        source(x, varArgs...) = [b...]
         f::Array{Float64,1} = -RapidFEM.assembleVector((source, initSoln), volAttrib, FeSpace,
        mesh, localReferenceSource!, problemDim, activeDimensions)
        #println("Fx = ", Fx)
        neumann(x; varArgs...) = begin
            center = [1.0, 0.5, 0.5]
            X,Y,Z = x-center
            r = sqrt(Y^2+Z^2)
            F_max = sqrt(Fy^2+Fz^2)
            F_max*[Fx, -Z/r,Y/r]
        end

        #println("initSoln = ", initSoln)
         f -= RapidFEM.assembleVector((neumann, initSoln), neumAttrib, FeSpace,
        mesh, localReferenceNeumann!, problemDim, activeDimensions)
        #println("external = ", f)

        hyperElasticParameters = (hyperModel, modelParams, initSoln)
         fσ::Array{Float64,1} = RapidFEM.assembleVector!(hyperElasticParameters, volAttrib, FeSpace,
        mesh, local_δE_S_Vector!, problemDim, activeDimensions)
        #println("\nnorm(fσ) = ", fσ,"\n")#, "\nnorm(initSoln) = ", initSoln)

         f+= fσ
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttrib1, mesh, problemDim)
        #println("f =", f)
        elasped = time() - start
        #println("Time for vector Assembly = ", elasped)
        return f
    end

    assemble_J(initSoln) = begin
        start = time()
        hyperElasticParameters = (hyperModel, modelParams, initSoln)
         J::SparseMatrixCSC = RapidFEM.assembleMatrix!(hyperElasticParameters, volAttrib, FeSpace,
        mesh, local_δE_Cᵀ_ΔE!, problemDim, activeDimensions)

         J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttrib1, mesh, problemDim)
        #println("J =", Matrix(J))
        elasped = time() - start
        #println("Time for Matrix Assembly = ", elasped)
        return J
    end

    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("cube2-Twist", mesh, [volAttrib], problemDim)
    #plot()
    for i ∈ 0:noOfSteps

        γ[i+1] = (i/noOfSteps)#^0.45

        Fx = γ[i+1]*Fx_max
        Fy = γ[i+1]*Fy_max
        Fz = γ[i+1]*Fz_max

        b .= γ[i+1]*b_max

        println("i = ", i, " Fx = ", Fx, " Fy = ", Fy, " Fz = ", Fz)

        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction1, dirchAttrib1, mesh, problemDim)

        #initSoln, convergenceData = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
        #    xtol = 1e-11, ftol = 1.e-5, relTol= 1e-8, iterations = 100, skipJacobian = 1 , printConvergence = true)
        solverResults = nlsolve(assemble_f, assemble_J, initSoln, show_trace = true, method = :trust_region, linesearch = BackTracking())
        initSoln .= solverResults.zero

        hyperElasticParameters = (hyperModel, modelParams, initSoln)

        SecondPiolaStress = RapidFEM.InvDistInterpolation([gaussianSecondPiolaStress],
        initSoln, [hyperElasticParameters],  FeSpace, mesh,  [volAttrib],
            problemDim, activeDimensions)

        DeformationGrad = RapidFEM.InvDistInterpolation([gaussianDeformationGrad],
         initSoln, [hyperElasticParameters],  FeSpace, mesh,  [volAttrib],
             problemDim, activeDimensions)

        GreenStrain = RapidFEM.InvDistInterpolation([gaussianGreenStrain],
        initSoln, [hyperElasticParameters],  FeSpace, mesh,  [volAttrib],
                problemDim, activeDimensions)

        CauchyStress = RapidFEM.InvDistInterpolation([gaussianCauchyStress],
        initSoln, [hyperElasticParameters],  FeSpace, mesh,  [volAttrib],
                        problemDim, activeDimensions)

        RapidFEM.vtkDataAdd!(vtkMeshData, (initSoln, SecondPiolaStress, GreenStrain, DeformationGrad, CauchyStress), (
        "Displacement", "Second Piola Stress", "Green Strain", "Deformation Gradient", "Cauchy Stress"), float(i), i)

        #plot(log10.(convergenceData.relNorm), linestyle= :dashdot)
    end
    #xlabel("Iterations")
    #ylabel("Log base 10 of Relative Convergence")
    RapidFEM.vtkSave(vtkMeshData)
end
