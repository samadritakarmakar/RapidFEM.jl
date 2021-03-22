using RapidFEM, LargeDeformations, SparseArrays, LinearAlgebra, PyPlot
include("hyperElasticLocalAssembly/hyperElastic.jl")

function hyperElasticity()
    #mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/HexahedralOrder2.msh")
    #mesh::Mesh = RapidFEM.readMesh("../test/MeshFiles/cube.msh")
    mesh::Mesh = RapidFEM.readMesh("../test/MeshFiles/cube1elmnt.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,4)
    neumAttrib::Tuple{Int64, Int64} = (2,3) #Force
    dirchAttrib1::Tuple{Int64, Int64} = (2,1) #Lock
    dirchAttrib2::Tuple{Int64, Int64} = (2,2) #move
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 10.0 #MPa
    ν::Float64 = 0.3
    #λ = (ν*E)/((1+ν)*(1-2*ν))
    #μ = E/(2*(1+ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))

    Fx::Float64 = 0.0
    Fy::Float64 = 0.0
    Fz::Float64 = 0.0

    Fx_max::Float64 = 0.1
    Fy_max::Float64 = 0.0
    Fz_max::Float64 = 0.0
    noOfSteps = 1

    b = [0.0, -0.5, 0.0]
    γ = zeros(noOfSteps)

    hyperModel = LargeDeformations.saintVenantModel
    #hyperModel = LargeDeformations.neoHookeanCompressibleModel
    modelParams::Tuple  = (λ, μ)

    totalDoF::Int64 = mesh.noOfNodes*problemDim
    initSoln::Array{Float64, 1} = zeros(totalDoF)

    DirichletFunction1(x; varArgs...) = zeros(problemDim)
    DirichletFunction2(x; varArgs...) = [0.0, (0.5+(x[2]-0.5)*cos(π/3)-(x[3]-0.5)*sin(π/3)-x[2])/2.0, (0.5+(x[2]-0.5)*sin(π/3)+(x[3]-0.5)*cos(π/3)-x[1])/2.0]
    assemble_f(initSoln) = begin
        start = time()
        source(x, varArgs...) = [b...]
         f::Array{Float64,1} = -RapidFEM.assembleVector((source, initSoln), volAttrib, FeSpace,
        mesh, localReferenceSource!, problemDim, activeDimensions)
        #println("Fx = ", Fx)
        neumann(x; varArgs...) = [Fx, Fy, Fz]

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
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttrib2, mesh, problemDim)
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
         J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttrib2, mesh, problemDim)
        #println("J =", Matrix(J))
        elasped = time() - start
        #println("Time for Matrix Assembly = ", elasped)
        return J
    end

    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("HyperElasticLagrange", mesh, [volAttrib], problemDim)
    plot()
    for i ∈ 1:noOfSteps

        γ[i] = (i/noOfSteps)#^0.45

        Fx = γ[i]*Fx_max
        Fy = γ[i]*Fy_max
        Fz = γ[i]*Fz_max

        println("i = ", i, " Fx = ", Fx, " Fy = ", Fy, " Fz = ", Fz)

        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction1, dirchAttrib1,
        mesh, problemDim)
        RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction2, dirchAttrib2,
        mesh, problemDim)

         initSoln, convergenceData = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
            xtol = 1e-11, ftol = 1.e-5, relTol= 1e-7, iterations = 100, skipJacobian = 1 , printConvergence = true)
        RapidFEM.vtkDataAdd!(vtkMeshData, (initSoln,), ("Displacement", ), float(i), i)
        plot(log10.(convergenceData.relNorm), linestyle= :dashdot)
    end
    xlabel("Iterations")
    ylabel("Natural Log of Relative Convergence")
    RapidFEM.vtkSave(vtkMeshData)
end
