using RapidFEM, SparseArrays, WriteVTK, FEMSparse

meshFile = "../../test/MeshFiles/cubeFullDirichletR3O1.msh"

#actualSolution(x) = return [sum(x.^2)+x[1]*x[2]*x[3]]
actualSolution(x) = [sum(sin.(x))]

function compareSolution(x)
    mesh::Mesh = RapidFEM.readMesh(meshFile)
    nodeKeys = collect(keys(mesh.Nodes))
    err = zeros(length(nodeKeys))
    for key ∈ nodeKeys
        err[key] = x[key] - actualSolution(mesh.Nodes[key])[1]
    end
    return err
end

function poissonEquation()
    #mesh::Mesh = RapidFEM.readMesh("../../test/OneElmntMsh/TetrahedralOrder2.msh")
    mesh::Mesh = RapidFEM.readMesh(meshFile)
    #mesh::Mesh = RapidFEM.readMesh("../../test/OneElmntMsh/HexahedralOrder1.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 1
    volAttrib::Tuple{Int64, Int64} = (3,2)
    dirchAttrib::Tuple{Int64, Int64} = (2,1)
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    parameterFunction(x, varArgs...) = [1.0]#, 1.0, 1.0]
    K::SparseMatrixCSC = RapidFEM.assembleMatrix(parameterFunction, volAttrib,
    FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)
    source(x, varArgs...) = [actualSolution(x)...]#, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib,
    FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    RapidFEM.applyDirichletBC!(f, K, actualSolution, dirchAttrib,
    mesh, problemDim)
    x::Vector = K\f
    err = compareSolution(x)
    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("poisson", mesh, [volAttrib],problemDim)
    RapidFEM.vtkDataAdd!(vtkMeshData, (x, err), ("Field","Error"))
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
end


