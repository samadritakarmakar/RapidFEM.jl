using RapidFEM, SparseArrays, WriteVTK

function poissonEquation()
    mesh::Mesh = RapidFEM.readMesh("../test/Bar.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,4)
    neumAttrib::Tuple{Int64, Int64} = (2,1)
    dirchAttrib::Tuple{Int64, Int64} = (2,2)
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    #println(mesh.Elements[3,3][1].nodeTags)
    K::SparseMatrixCSC = RapidFEM.assembleMatrix(volAttrib, FeSpace, mesh, RapidFEM.local_lagrange_K, problemDim, activeDimensions)
    source(x) = [0.0, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib, FeSpace, mesh, RapidFEM.localSource, problemDim, activeDimensions)
    neumann(x) = [0.0, 0.1, 0.0]
    f += RapidFEM.assembleVector(neumann, neumAttrib, FeSpace, mesh, RapidFEM.localNeumann, problemDim, activeDimensions)
    DirichletFunction(x) = zeros(problemDim)
    RapidFEM.applyDirichletBC!(K, f, DirichletFunction, dirchAttrib, mesh, problemDim)
    x::Vector = K\f
    vtkfile = RapidFEM.InitializeVTK(x, "field",mesh, [volAttrib], problemDim)
    vtkfile["Displacement"] = x
    RapidFEM.vtkSave(vtkfile)
end
