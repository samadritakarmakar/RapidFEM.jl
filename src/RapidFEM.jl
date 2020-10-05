__precompile__()
module RapidFEM
using FEMSparse, SparseArrays, LinearAlgebra, WriteVTK
include("FEM/elements.jl")
include("ShapeFunctions/shapeFunction.jl")
include("FEM/boundaryCondition.jl")
include("FEM/assembly.jl")
include("FEM/feSpace.jl")
include("FEM/dofUtils.jl")
include("FEM/postProcess.jl")
include("FEM/utils.jl")
include("Models/general.jl")
include("Models/linearElasticity.jl")
include("Models/convectionFluid.jl")
include("Output/WriteToVTK.jl")


#From FEM
    export AbstractElement, LineElement, TriElement, QuadElement, TetElement,HexElement, shapeFunction, ipPoint
#feSpace
    export createFeSpace, feSpace!, lagrange, get_∂x_∂ξ, getFunction_dΩ, getFunction_dS, getFunction_dL, getFunction_∂ξ_∂x, getInterpolated_x, getNodes, getVectorNodes
#boundaryCondition
    export applyDirichletBC!, assembleVector, assembleMatrix
#utils
    export elmntSizeAlongVel
#postProcess
    export InvDistInterpolation, voigtToTensor
#From Mesh
export Mesh, readMesh, getNoOfElements, getCoordArray
#From Quadrature
export getQuadrature
#From ShapeFunction
export IpPoint, ShapeFunction, calculateShapeFunctions
#From Models
#general
    export local_∇v_λ_∇u, local_v_ρ_u, localSource, localNeumann, localScalar, localScalarNeumann
#linearElasticity
    export local_∇v_C_∇u, createVoigtElasticTensor, getTensorMapping, gaussianStress
#convectionFluid
    export local_v_λ_∇u_Vector, local_v_λ_∇u_Scalar
#From WriteToVTK
    export VTKMeshData, InitializeVTK, vtkSave, vtkDataAdd
end # module
