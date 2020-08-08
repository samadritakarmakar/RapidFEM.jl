__precompile__()
module RapidFEM
using FEMSparse, SparseArrays, LinearAlgebra, WriteVTK
include("FEM/elements.jl")
include("FEM/feSpace.jl")
include("FEM/dofUtils.jl")
include("Models/general.jl")
include("Models/linearElasticity.jl")
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
#general
export local_lagrange_K, localSource, localNeumann, localScalar, localScalarNeumann
#linearElasticity
export local_∇v_C_∇u, createVoigtElasticTensor, getTensorMapping
#From WriteToVTK
export InitializeVTK, vtkSave

end # module
