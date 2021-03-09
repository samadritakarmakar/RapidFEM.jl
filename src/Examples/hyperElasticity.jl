using RapidFEM, LargeDeformations, SparseArrays
include("hyperElasticLocalAssembly/hyperElastic.jl")

function hyperElasticity()
    mesh::Mesh = RapidFEM.readMesh("../test/OneElmntMsh/HexahedralOrder1.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,2) #Force
    dirchAttrib::Tuple{Int64, Int64} = (2,1) #Lock
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    Fx::Float64 = 0.5

    hyperModel = LargeDeformations.saintVenantModel
    modelParams::Tuple  = (λ, μ)

    totalDoF::Int64 = mesh.noOfNodes*problemDim
    initSoln::Array{Float64, 1} = zeros(totalDoF)

    DirichletFunction(x; varArgs...) = zeros(problemDim)

    assemble_f(initSoln) = begin
        source(x, varArgs...) = [0.0, 0.0, 0.0]
        f::Array{Float64,1} = -RapidFEM.assembleVector((source, initSoln), volAttrib, FeSpace,
        mesh, localCurrentSource!, problemDim, activeDimensions)
        println("Fx = ", Fx)
        neumann(x; varArgs...) = [Fx, 0.0, 0.0]

        f -= RapidFEM.assembleVector((neumann, initSoln), neumAttrib, FeSpace,
        mesh, localReferenceNeumann!, problemDim, activeDimensions)
        #println("external = ", f)

        hyperElasticParameters = (hyperModel, modelParams, initSoln)
        fσ::Array{Float64,1} = RapidFEM.assembleVector!(hyperElasticParameters, volAttrib, FeSpace,
        mesh, local_∇v_σ_Vector!, problemDim, activeDimensions)
        RapidFEM.applyNLDirichletBC_on_f!(fσ, dirchAttrib, mesh, problemDim)
        #println("\nnorm(fσ) = ", fσ,"\n")#, "\nnorm(initSoln) = ", initSoln)

        f+= fσ
        RapidFEM.applyNLDirichletBC_on_f!(f, dirchAttrib, mesh, problemDim)
        println("f =", f)
        return f
    end

    assemble_J(initSoln) = begin
        hyperElasticParameters = (hyperModel, modelParams, initSoln)
        J::SparseMatrixCSC = RapidFEM.assembleMatrix!(hyperElasticParameters, volAttrib, FeSpace,
        mesh, local_∇v_Cᵀ_∇u!, problemDim, activeDimensions)
        J = RapidFEM.applyNLDirichletBC_on_J!(J, dirchAttrib, mesh, problemDim)
        return J
    end

    RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttrib,
    mesh, problemDim)

    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("HyperElastic", mesh, [volAttrib], problemDim)
    RapidFEM.applyNLDirichletBC_on_Soln!(initSoln, DirichletFunction, dirchAttrib,
    mesh, problemDim)
    initSoln = RapidFEM.simpleNLsolve(assemble_f, assemble_J, initSoln;
        xtol = 1e-12, ftol = 1e-8, iterations = 2000, skipJacobian =1 , printConvergence = true)
    RapidFEM.vtkDataAdd!(vtkMeshData, (initSoln,), ("Displacement", ), float(i), i)

    RapidFEM.vtkSave(vtkMeshData)
end
