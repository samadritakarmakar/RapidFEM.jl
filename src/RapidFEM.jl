module RapidFEM
using FEMSparse, SparseArrays, LinearAlgebra
include("FEM/elements.jl")
include("FEM/feSpace.jl")
include("FEM/dofUtils.jl")
include("Models/general.jl")


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
end # module
