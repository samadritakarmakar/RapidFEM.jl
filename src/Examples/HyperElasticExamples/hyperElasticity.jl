using RapidFEM, LargeDefs, SparseArrays, Tensors, PyPlot
include("../hyperElasticLocalAssembly/hyperElastic.jl")

function hyperElasticity()
    #mesh::Mesh = RapidFEM.readMesh("../../test/OneElmntMsh/HexahedralOrder1.msh")
    mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/cubeTest.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,5)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 10.0 #MPa
    ν::Float64 = 0.3
    #λ = (ν*E)/((1+ν)*(1-2*ν))
    #μ = E/(2*(1+ν))
    μ = 3.8
    λ = 5.0

    println("λ = ", λ, " μ = ", μ)

    Fx::Float64 = 2.5

    #hyperModel = LargeDefs.neoHookeanCompressible
    hyperModel = LargeDefs.neoHookean
    modelParams::Tuple  = (λ, μ)

    totalDoF::Int64 = mesh.noOfNodes*problemDim
    initSoln::Array{Float64, 1} = zeros(totalDoF)

    DirichletFunction(x; varArgs...) = zeros(problemDim)

    assemble_f(initSoln) = begin
        source(x, varArgs...) = [0.0, 0.0, 0.0]
        f::Array{Float64,1} = -RapidFEM.assembleVector((source, initSoln), volAttrib, FeSpace,
        mesh, localReferenceSource!, problemDim, activeDimensions)
        neumann(x; varArgs...) = [Fx, 0.0, 0.0]

        f -= RapidFEM.assembleVector((neumann, initSoln), neumAttrib, FeSpace,
        mesh, localReferenceNeumann!, problemDim, activeDimensions)
        #println("external = ", f)

        hyperElasticParameters = (hyperModel, modelParams, initSoln)
         fσ::Array{Float64,1} = RapidFEM.assembleVector!(hyperElasticParameters, volAttrib, FeSpace,
        mesh, local_δE_S_Vector!, problemDim, activeDimensions)
        #RapidFEM.applyNLDirichletBC_on_f!(fσ, dirchAttrib, mesh, problemDim)
        #println("\nnorm(fσ) = ", fσ,"\n")#, "\nnorm(initSoln) = ", initSoln)

        f+= fσ
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttrib, mesh, problemDim)
        #println("f =", f)
        return f
    end

    assemble_J(initSoln) = begin
        hyperElasticParameters = (hyperModel, modelParams, initSoln)
        J::SparseMatrixCSC = RapidFEM.assembleMatrix!(hyperElasticParameters, volAttrib, FeSpace,
        mesh, local_δE_Cᵀ_ΔE!, problemDim, activeDimensions)
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttrib, mesh, problemDim)
        return J
    end

    i = 1

    RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttrib,
    mesh, problemDim)
    solnNode = 3
    
    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("HyperElastic", mesh, [volAttrib], problemDim)
    RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttrib,
    mesh, problemDim)
    initSoln, convergenceData = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
       xtol = 1e-11, ftol = 1.e-5, relTol= 1e-8, iterations = 100, skipJacobian = 1 , printConvergence = true)
       println("Soln at node $solnNode = ", initSoln[problemDim*(solnNode-1)+1:problemDim*solnNode])
    RapidFEM.vtkDataAdd!(vtkMeshData, (initSoln,), ("Displacement", ), float(i), i)
    plot(log10.(convergenceData.relNorm), linestyle= :dashdot)

    xlabel("Iterations")
    ylabel("Log base 10 of Relative Convergence")
    RapidFEM.vtkSave(vtkMeshData)
end
