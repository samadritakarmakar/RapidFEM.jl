#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

__precompile__()
module RapidFEM
using FEMSparse, FastGaussQuadrature, SparseArrays, LinearAlgebra, Arpack, WriteVTK
include("FEM/elements.jl")
include("Mesh/mesh.jl")
include("Mesh/aspectRatio.jl")
include("Mesh/meshExtra.jl")
include("Mesh/addBoundaryElements.jl")
include("Mesh/extrudeSurfaceToVolume.jl")
include("Mesh/updateMesh.jl")
include("Quadrature/quadrature.jl")
include("ShapeFunctions/langrange.jl")
include("ShapeFunctions/shapeFunction.jl")
include("FEM/boundaryCondition.jl")
include("FEM/freeVibration.jl")
include("FEM/assembly.jl")
include("FEM/feSpace.jl")
include("FEM/coupledComponents.jl")
include("FEM/dofUtils.jl")
include("FEM/postProcess.jl")
include("FEM/SprLikeRecovery.jl")
include("FEM/utils.jl")
include("FEM/normals.jl")
include("NonLinear/nonLinearUtils.jl")
include("NonLinear/simpleNLsolve.jl")
include("Dynamic/singleStep_pj.jl")
include("Dynamic/newmark.jl")
include("Models/general.jl")
include("Models/linearElasticity.jl")
include("Models/convectionFluid.jl")
#include("Models/smallStrainPlasticity.jl")
include("Output/SolutionVectorIn3d.jl")
include("Output/WriteToVTK.jl")
include("Output/WriteToMesh.jl")



#From FEM
    export AbstractElement, LineElement, TriElement, QuadElement, TetElement, HexElement, PointElement, shapeFunction, ipPoint
##feSpace
    export createFeSpace, feSpace!, lagrange, get_∂x_∂ξ, getFunction_dΩ, getFunction_dS, getFunction_dL, getFunction_∂ξ_∂x, getInterpolated_x
    export get_dΩ, get_∂ξ_∂x, get_dS, get_dL, get_∂ϕ_∂x, get_ϕ, getNoOfElementNodes, getNoOfElementIpPoints
#dofUtils
    export getNodes, getVectorNodes, getSolAtElement, getCoupledGlobalSols, getMixedGlobalSols
#CoupledComponents
    export CoupledComponents
##boundaryCondition
    export applyDirichletBC!, applyDynamicDirichletBC!
    export applyNLDirichletBC_on_J!, applyNLDirichletBC_on_Soln!
    export applyNLDirichletBC_on_f!
    export getUniqueNodes
##freeVibration
    export applyFreeVibrationBC!, solveFreeVibration, getFreeVibrationDisplacement, solveFreeVibrationDense, scaleModeShapes!
##assembly
    export assembleVector, assembleMatrix, assembleScalar
    export assembleVector!, assembleMatrix!, assembleScalar!
##utils
    export elmntSizeAlongVel, get_∂u_∂x!, get_∂u_∂x, get_u!, get_u, getCurrentCoordArray!, getCurrentCoordArray
##postProcess
    export InvDistInterpolation, voigtToTensor
##SprLikeRecovery
    export SprLikeRecovery
##normals
    export getElementSurfaceNormal, getSurfaceNormals
#From Mesh
    #mesh.jl
    export Mesh, readMesh, getNoOfElements, getCoordArray, updateNodePositions!, getCurrentNodeDict
    export getElementTypeAndOrder, createNewElement, getNodesDict, getTotalElements
    #aspectRatio.jl
    export getAspectRatioOfElement, getAspectRatios, getAspectRatioOfTriElement, getAspectRatioOfTetElement, getAspectRatioOfQuadElement
    #meshExtra
    export MeshExtra, getAttributesInDimension, getNodeToElementMap
    export getNodeToElementMap!, getAllFaces, getBoundaryFaces
    export getAllMaterialBoundaryFaces, getAllBoundaryNodes, getAllInternalNodes
    #addBoundaryElements
    export addBoundaryElements!, removeElements!
    #updateMesh
    export replaceAndAddElements!
    #extrudeSurfaceToVolume
    export tapereExtrusionFunction, extrudeMeshSurfaceToVolume
#From Quadrature
export gauss
#From ShapeFunction
export IpPoint, ShapeFunction, calculateShapeFunctions
#From NonLinear
##simpleNLsolve
export simpleNLsolve, Convergence
##nonLinearUtils
export berganIncrement
#From Dynamic
##singleStep_pj
export get_SSpj_A_meanU_f, SSpj_getFinal_A_b, updateSolution!, update_f!
##newmark
export getNewmark_Disp_Velocity, ImplicitNewmark, getNewmarkAccsFromDispVel
##Generalized Alpha
export getGenAlpha_DispVelAcc, getImplicitGenAlphaParameters, getExplicitGenAlphaParameters
export getGenAlphaParameters, GeneralizedAlpha
export getHhtAlphaParameters, getHhtAlpha_DispVecAcc
#From Models
##general
    export local_∇v_λ_∇u!, local_v_ρ_u!, localBoundary_v_ρ_u!, localSource!, localNeumann!, localScalar!, localScalarNeumann!
##linearElasticity
    export local_∇v_C_∇u!, createVoigtElasticTensor, getTensorMapping, gaussianStress, getVoigtIndex
##convectionFluid
    export local_v_λ_∇u_Vector!, local_v_λ_∇u_Scalar!
##SmallStrainPlastic
#    export local_∇v_Cᵀ_∇u!, local_∇v_σ_Vector!, j2Model, initParams_j2
#    export updateStateDict4rmBuffer
#From SolutionVectorIn3d
    export getSolutionVectorIn3d
#From WriteToVTK
    export VTKMeshData, InitializeVTK, vtkSave, vtkDataAdd!, createVTKFile
#From WriteToMesh
    export getGmshElementTypeNo, writeMesh
end # module
