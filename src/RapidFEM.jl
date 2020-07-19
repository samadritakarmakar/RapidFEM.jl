module RapidFEM
include("FEM/elements.jl")
include("FEM/feSpace.jl")
#From FEM
export AbstractElement, LineElement, TriElement, QuadElement, TetElement, HexElement, shapeFunction, ipPoint
export createFeSpace, feSpace!, lagrange, get_∂x_∂ξ, getFunction_dΩ, getFunction_∂ξ_∂x, getInterpolated_x
#From Mesh
export Mesh, readMesh, getNoOfElements, getCoordArray
#From Quadrature
export getQuadrature
#From ShapeFunction
export IpPoint, ShapeFunction, calculateShapeFunctions
end # module
