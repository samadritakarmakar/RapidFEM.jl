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
    AttributeName::Dict{Tuple{Int64, Int64},String}
    noOfNodes::Int64
    Nodes::Dict{Int64, Array{Float64, 1}}
    noOfElements::Int64
    Elements::Dict{Tuple{Int64, Int64}, Array{AbstractElement, 1}}
    meshSoftware::String
end

"""The readMesh function is responsible to read the mesh and store the data in
a mesh struct type. It has been designed to be versatile to read different kinds
mesh files but as of now it only support Gmsh ASCII files of version 2.02

    mesh::Mesh = readMesh("../test/Bar.msh")
"""
function readMesh(meshFileName::String)
    attributes::Array{Tuple{Int64, Int64},1} = []
    AttributeName::Dict{Tuple{Int64, Int64},String} = Dict{Tuple{Int64, Int64},String}()
    Nodes::Dict{Int64, Array{Float64, 1}} = Dict{Int64, Array{Float64, 1}}()
    Elements::Dict{Tuple{Int64, Int64}, Array{AbstractElement, 1}}=
     Dict{Tuple{Int64, Int64}, Array{AbstractElement, 1}}()
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

"""Updates Nodal Positions, given that the change in it's positions are known.

    updateNodePositions!(mesh::Mesh, changeInPosition::Array{Float64, 1}, activeDimensions = [1,1,1])
"""
function updateNodePositions!(mesh::Mesh, changeInPosition::Array{Float64, 1}, activeDimensions = [1,1,1])
    nodes = sort(collect(keys(mesh.Nodes)))
    rangeDim = createDimRange()
    problemDim = sum(activeDimensions)
    range = getRange(rangeDim, activeDimensions)
    for node ∈ nodes
        mesh.Nodes[node][range] = mesh.Nodes[node][range] + changeInPosition[problemDim*(node-1)+1:problemDim*node]
    end
    return nothing
end

"""Returns the element type and order from the dimension of the element and the length of its node tags.

        elementType, order = getElementTypeAndOrder(dimension, newElementNodeTags)
"""
function getElementTypeAndOrder(dimension::Int64, newElementNodeTags::Vector{Int64})
    lengthElementNodeTags = length(newElementNodeTags)
    elementType = TriElement
    order = 1
    notFound = false
    if dimension == 0
        elementType, order = PointElement, 1
    elseif dimension == 1
        elementType = LineElement
        if lengthElementNodeTags == 2
            order = 1
        elseif lengthElementNodeTags == 3
            order = 2
        elseif lengthElementNodeTags == 4
            order = 3
        else
            notFound = true
        end
    elseif dimension == 2
        if lengthElementNodeTags == 3 
            elementType, order = TriElement, 1
        elseif lengthElementNodeTags == 6
            elementType, order = TriElement, 2
        elseif lengthElementNodeTags == 10
            elementType, order = TriElement, 3
        elseif lengthElementNodeTags == 4
            elementType, order = QuadElement, 1
        elseif lengthElementNodeTags == 9
            elementType, order = QuadElement, 2
        elseif lengthElementNodeTags == 16
            elementType, order = QuadElement, 3
        else
            notFound = true
        end
    elseif dimension == 3 
        if lengthElementNodeTags == 4 
            elementType, order = TetElement, 1
        elseif lengthElementNodeTags == 10
            elementType, order = TetElement, 2
        elseif lengthElementNodeTags == 20
            elementType, order = TetElement, 3
        elseif lengthElementNodeTags == 8
            elementType, order = HexElement, 1
        elseif lengthElementNodeTags == 27
            elementType, order = HexElement, 2
        elseif lengthElementNodeTags == 64
            elementType, order = HexElement, 3
        else
            notFound = true
        end
    else
        notFound = true           
    end
    if notFound
        error("Element with attribute $dimension and node length $lengthElementNodeTags not supported.")
    end
    return elementType, order
end

"""Returns a new Element from with a given label and attribute. Recognizes type of element and order from dimension in attib[1] and 
length of newElementNodeTags

    createNewElement(attrib, label, newElementNodeTags)
"""
function createNewElement(attrib::Tuple{Int64, Int64}, label::Int64, newElementNodeTags::Vector{Int64})
    
    lengthElementNodeTags = length(newElementNodeTags)
    elementType, order = getElementTypeAndOrder(attrib[1], newElementNodeTags)
    return elementType(label, [attrib[2], attrib[2]], 
        newElementNodeTags, lengthElementNodeTags, order)   
end

function getNodesDict(primaryNodeArray::AbstractMatrix{Float64}, activeDims::Vector{Int64})
    Nodes = Dict{Int64, Array{Float64, 1}}()
    coord = zeros(3)
    #for nodeNo ∈ 1:size(primaryNodeArray,2)
    nodeNo = 1
    for nodeCoord ∈ eachcol(primaryNodeArray)
       activeDimNo = 1
       dimNo = 1
       for activeDim ∈ activeDims
          coord[dimNo] = activeDim == 1 ? nodeCoord[activeDimNo] : 0.0
          activeDimNo += activeDim == 1 ? 1 : 0
          dimNo +=1
       end
       Nodes[nodeNo] = deepcopy(coord)
       nodeNo += 1
    end
    return Nodes
 end

 function getTotalElements(Elements::Dict{Tuple{Int64, Int64}, Array{AbstractElement, 1}})
    totalElements = 0
    for attrib ∈ keys(Elements)
        totalElements += length(Elements[attrib])
    end
    return totalElements
end

function getTotalElements(attribute::Tuple{Int64, Int64}, mesh::Mesh)
    if attribute ∈ keys(mesh.Elements)
        return length(mesh.Elements[attribute])
    end
    return 0
end

function getTotalElements(mesh::Mesh)
    totalElements = 0
    for attribute ∈ mesh.attributes
        totalElements += getTotalElements(attribute, mesh)
    end
    return totalElements
end
