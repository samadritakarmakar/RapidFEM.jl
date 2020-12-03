#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
 
include("gmshReader.jl")
"""Mesh struct keeps the data regarding the mesh and
it's properties."""
mutable struct Mesh
    noOfAttrib::Int64
    attributes::Array{Tuple{Int64, Int64},1}
    AttributeName::Dict{Any, Any}
    noOfNodes::Int64
    Nodes::Dict{Any,Any}
    noOfElements::Int64
    Elements::Dict{Any,Any}
    meshSoftware::String
end

"""The readMesh function is responsible to read the mesh and store the data in
a mesh struct type. It has been designed to be versatile to read different kinds
mesh files but as of now it only support Gmsh ASCII files of version 2.02

    mesh::Mesh = readMesh("../test/Bar.msh")
"""
function readMesh(meshFileName::String)
    attributes::Array{Tuple{Int64, Int64},1} = []
    AttributeName::Dict{Any,Any} = Dict()
    Nodes::Dict{Any,Any} = Dict()
    Elements::Dict{Any, Any} = Dict()
    noOfAttrib::Int64 = 0
    noOfNodes::Int64 = 0
    noOfElements::Int64 = 0
    meshSoftware::String = ""
    posExnt = findlast('.', meshFileName)
    if meshFileName[posExnt:end] == ".msh"
        noOfAttrib, noOfNodes, noOfElements = readGmshFile!(meshFileName, attributes,
            AttributeName, Nodes, Elements)
            meshSoftware = "gmsh"
    else
        error("This mesh type is not yet supported.")
    end
    mesh::Mesh = Mesh(noOfAttrib, attributes, AttributeName, noOfNodes, Nodes, noOfElements, Elements, meshSoftware)
    return mesh
end

"""This function returns the number of elements of a certain set of elements
for a given atttribute.

    noOfElements::Int64 = getNoOfElements(mesh, attribute)
"""
function getNoOfElements( mesh::Mesh, attribute::Tuple{Int64, Int64})
    return length(mesh.Elements[attribute])
end

"""This function extracts the Coordinate array for a certain element.
As an example, the output for a 1st order element would look like the
following:

[x₁, x₂, x₃;
y₁, y₂, y₃;
z₁, z₂, z₃;]

    coordArray::Array{Float64,2} = getCoordArray(mesh, element)
"""
function getCoordArray(mesh::Mesh,element::AbstractElement)::Array{Float64,2}
    #nodeTags::Array{Int64} = mesh.Elements[attribute][elementNo].nodeTags
    #noOfElementNodes::Int64 = mesh.Elements[attribute][elementNo].noOfElementNodes
    nodeTags::Array{Int64} = element.nodeTags
    noOfElementNodes::Int64 = element.noOfElementNodes
    CoordArray::Array{Float64,2} = Array{Float64}(undef, 3, noOfElementNodes)
    for elmntNodeNum::Int64 ∈ 1:noOfElementNodes
        CoordArray[:,elmntNodeNum] = mesh.Nodes[nodeTags[elmntNodeNum]]
    end
    return CoordArray
end
