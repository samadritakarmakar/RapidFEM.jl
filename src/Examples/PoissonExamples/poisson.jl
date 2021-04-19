#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
using RapidFEM, SparseArrays, WriteVTK, FEMSparse

function poissonEquation()

    #This function is used to load the mesh File.

    #mesh::Mesh = RapidFEM.readMesh("../../test/OneElmntMsh/TetrahedralOrder2.msh")
    mesh::Mesh = RapidFEM.readMesh("../../test/MeshFiles/BarHex.msh")
    #mesh::Mesh = RapidFEM.readMesh("../../test/OneElmntMsh/HexahedralOrder1.msh")
    
    # This is used to store the Shape Functions, it's 
    #gradients and hessinans at every Integration point 
    #for a certain type of element and it's order.
    FeSpace = RapidFEM.createFeSpace()
    
    #This represents the degrees of Freedom of the solution in each Node.
    problemDim::Int64 = 1

    #The Volume attribute represents the the attribute used to identify the attribute of the domain
    #being integrated over. The first entry is of the dimension of the whole domain,
    #the second is of the domain atttribute.
    volAttrib::Tuple{Int64, Int64} = (3,3)
    #The neumann attribute represents the the attribute used to identify the attribute of the surface
    #being integrated over. The first entry is of the dimension of the surface domain,
    #the second is of the surface atttribute.
    neumAttrib::Tuple{Int64, Int64} = (2,2)
    #The dirichlet attribute represents the the attribute used to identify the attribute of the surface
    #being set as dirichlet boundary.  The first entry is of the dimension of the 
    #surface domain, the second is of the surface atttribute.
    dirchAttrib::Tuple{Int64, Int64} = (2,1)
    #Tells the Library of the relevant dimensions in the mesh that are to be considered.
    #[1,1,1] says that x,y & z axis are the dimensions to be considered.
    #[1, 1, 0] says that x,y axis are the dimensions to be considered.
    #[0, 1, 1] says that y,z axis are the dimensions to be considered.
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    #This function is sent to the Matrix assembler used to assemble the stiffness or 
    #mass matrix. This function can be used to fix coeffients, such as conductivity values
    #for heat transfer problems. 
    parameterFunction(x, varArgs...) = [1.0]#, 1.0, 1.0]
    #This function returns the stiffness or mass matrix. It accepts arguments such as the 
    #parameter function the volume attribute, the created FeSpace, the mesh, the local
    #matrix assembler and the active dimensions. The local matrix assembler can be custom made. 
    K::SparseMatrixCSC = RapidFEM.assembleMatrix(parameterFunction, volAttrib,
    FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)
    source(x, varArgs...) = [0.0]#, 0.0, 0.0]
    
    f::Vector = RapidFEM.assembleVector(source, volAttrib,
    FeSpace, mesh, RapidFEM.localSource!, problemDim, activeDimensions)
    neumann(x, varArgs...) = [0.1]#, 0.0, 0.0]
    f += RapidFEM.assembleVector(neumann, neumAttrib,
    FeSpace, mesh, RapidFEM.localNeumann!, problemDim, activeDimensions)
    DirichletFunction(x, varArgs...) = zeros(problemDim)
    K = RapidFEM.applyDirichletBC!(f, K, DirichletFunction, dirchAttrib,
    mesh, problemDim)
    x::Vector = K\f
    vtkMeshData::VTKMeshData = RapidFEM.InitializeVTK("poisson", mesh, [volAttrib],problemDim)
    RapidFEM.vtkDataAdd!(vtkMeshData, (x,), ("Field",))
    RapidFEM.vtkSave(vtkMeshData)
    return nothing
end
