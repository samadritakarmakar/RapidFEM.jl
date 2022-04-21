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

function getElementTypeAndOrder(attrib::Tuple{Int64, Int64}, newElementNodeTags::Vector{Int64})
    lengthElementNodeTags = length(newElementNodeTags)
    elementType = TriElement
    order = 1
    notFound = false
    if attrib[1] == 0
        elementType, order = PointElement, 1
    elseif attrib[1] == 1
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
    elseif attrib[1] == 2
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
    elseif attrib[1] == 3 
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
        error("Element with attribute $attrib and node length $lengthElementNodeTags not supported.")
    end
    return elementType, order
end

"""Returns a new Element from with a given label and attribute. Recognizes type of element and order from dimension in attib[1] and 
length of newElementNodeTags

    createNewElement(attrib, label, newElementNodeTags)
"""
function createNewElement(attrib::Tuple{Int64, Int64}, label::Int64, newElementNodeTags::Vector{Int64})
    
    lengthElementNodeTags = length(newElementNodeTags)
    elementType, order = getElementTypeAndOrder(attrib, newElementNodeTags)
    return elementType(label, [attrib[2], attrib[2]], 
        newElementTags, lengthElementNodeTags, order)   
end

"""Replaces and Adds Elements in a Set of element attibutes. If badElAttribsElNos is used directly then it is recommended 
that all element attributes be the same.

    replaceAndAddElements!(mesh, badElAttribsElNos, newElementNodeTags)

    replaceAndAddElements!(mesh, attrib, elementNos, newElementNodeTags)
"""

function replaceAndAddElements!(mesh::Mesh, 
    badElAttribsElNos::Set{Tuple{Tuple{Int64, Int64}, Int64}},
    newElementNodeTags::Union{Vector{Vector{Int64}}, Set{Vector{Int64}}})

    if newElementNodeTags isa Set
        newElementNodeTags = collect(newElementNodeTags)
    end

    lengthNewElementTags = length(newElementNodeTags)
    newElementTagNo = 1
    deletedElements = 0
    attrib = (0,0)
    for badElAttribsElNo ∈ badElAttribsElNos
        attrib, changeTagElNo = badElAttribsElNo
        if newElementTagNo <= lengthNewElementTags
            label = mesh.Elements[attrib][changeTagElNo].label
            element = createNewElement(attrib, label, newElementNodeTags[newElementTagNo])
            
            mesh.Elements[attrib][changeTagElNo] = element
            newElementTagNo += 1
        else
            deleteat!(mesh.Elements[attrib],changeTagElNo)
            deletedElements += 1
        end
    end
    addElements = 0
    while newElementTagNo <= lengthNewElementTags
        newLabel = mesh.Elements[attrib][end].label + 1 
        element = createNewElement(attrib, newLabel, newElementNodeTags[newElementTagNo])
            push!(mesh.Elements[attrib], element)
       
        newElementTagNo += 1
        addElements += 1
    end
    mesh.noOfElements += addElements - deletedElements
end

function replaceAndAddElements!(mesh::Mesh, attrib::Tuple{Int64, Int64},
    elementNos::Union{Vector{Int64}, Set{Int64}}, 
    newElementNodeTags::Union{Vector{Vector{Int64}}, Set{Vector{Int64}}})

    badElAttribsElNos = Set{Tuple{Tuple{Int64, Int64}, Int64}}()
    for elementNo ∈ elementNos
        push!(badElAttribsElNos, (attrib, elementNo))
    end
    replaceAndAddElements!(mesh, badElAttribsElNos, newElementNodeTags)
end