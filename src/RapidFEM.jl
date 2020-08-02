__precompile__()
module RapidFEM
using FEMSparse, SparseArrays, LinearAlgebra, WriteVTK
include("FEM/elements.jl")
include("FEM/feSpace.jl")
include("FEM/dofUtils.jl")
include("Models/general.jl")
include("Output/WriteToVTK.jl")


#From FEM
export AbstractElement, LineElement, TriElement, QuadElement, TetElement,HexElement, shapeFunction, ipPoint
export createFeSpace, feSpace!, lagrange, get_∂x_∂ξ, getFunction_dΩ, getFunction_dS, getFunction_dL, getFunction_∂ξ_∂x, getInterpolated_x, getNodes, getVectorNodes
export applyDirichletBC!, assembleVector, assembleMatrix
#From Mesh
export Mesh, readMesh, getNoOfElements, getCoordArray
#From Quadrature
export getQuadrature
#From ShapeFunction
export IpPoint, ShapeFunction, calculateShapeFunctions
#From Models
export local_lagrange_K, localSource, localNeumann
#From WriteToVTK
export InitializeVTK, vtkSave

###################Precompiling Functions#################
function precompileEquations()
    mesh::Mesh = RapidFEM.readMesh("src/test/Bar.msh")
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
    vtkfile = RapidFEM.InitializeVTK(x, "/tmp/field", mesh, [volAttrib], problemDim)
    vtkfile["field"] = x
    RapidFEM.vtkSave(vtkfile)
end
end # module
