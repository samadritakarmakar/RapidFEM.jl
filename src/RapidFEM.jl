module RapidFEM
include("FEM/elements.jl")
include("FEM/feSpace.jl")
include("FEM/dofUtils.jl")

#From FEM
export AbstractElement, LineElement, TriElement, QuadElement, TetElement, HexElement, shapeFunction, ipPoint
export createFeSpace, feSpace!, lagrange, get_∂x_∂ξ, getFunction_dΩ, getFunction_dS, getFunction_dL, getFunction_∂ξ_∂x, getInterpolated_x, getNodes, getVectorNodes
#From Mesh
export Mesh, readMesh, getNoOfElements, getCoordArray
#From Quadrature
export getQuadrature
#From ShapeFunction
export IpPoint, ShapeFunction, calculateShapeFunctions
end # module
