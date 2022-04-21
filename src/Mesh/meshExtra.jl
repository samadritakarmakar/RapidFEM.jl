using Combinatorics

mutable struct MeshExtra
    nodeToElementMap::Dict{Int64, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}
    allFaces::Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}
    #boundaryFaces::Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}
    boundaryNodes::Array{Int64, 1}
    internalNodes::Array{Int64, 1}
    elementNodeTagsToElementMap::Dict{Vector{Int64}, Tuple{Tuple{Int64, Int64}, Int64}}
    #MeshExtra(nodeToElementMap::Dict{Int64, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}) = new(nodeToElementMap)
end

function MeshExtra()
    MeshExtra(Dict{Int64, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}(), 
    Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}(),
    zeros(Int64, 0), zeros(Int64, 0), 
    Dict{Vector{Int64}, Tuple{Tuple{Int64, Int64}, Int64}}())
end

function MeshExtra(mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1})
    nodeToElementMap = getNodeToElementMap(mesh, attributeArray)
    allFaces = getAllFaces(mesh, attributeArray)
    boundaryFaces = getBoundaryFaces(allFaces)
    boundaryNodes = getAllBoundaryNodes(boundaryFaces)
    internalNodes = getAllInternalNodes(boundaryNodes, mesh)
    elNodeTagsToElMap = getElNodeTagsToElementMap(mesh::Mesh, attributeArray)
    MeshExtra(nodeToElementMap, allFaces, boundaryNodes, internalNodes, elNodeTagsToElMap)
end

"""Retrives all attibutes belonging to a certain dimension"""
function getAttributesInDimension(mesh::Mesh, dimension::Int64)
    attributeKeys = keys(mesh.Elements)
    dimAttribsArray = Tuple{Int64,Int64}[]
    for attributeKey ∈ attributeKeys
        if attributeKey[1] == dimension
            push!(dimAttribsArray, attributeKey)
        end
    end
    return dimAttribsArray
end

"""Computes and creates extra mesh data from mesh and it's given dimension"""
function MeshExtra(mesh::Mesh, meshDimension::Int64 = 3)
    attributeArray = getAttributesInDimension(mesh, meshDimension)
    return MeshExtra(mesh, attributeArray)
end

function getNodeToElementMap(mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1})
    
    nodeToElementMap = Dict{Int64, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}()
    for attribute ∈ attributeArray
        elements = mesh.Elements[attribute...]
        noOfElements = size(elements,1)
        for elementNo ∈ 1:noOfElements
            element = elements[elementNo]
            for nodeTag ∈ element.nodeTags
            #for nodeNo ∈ 1:length(element)
                #nodeTag = element.nodeTags[nodeNo]
                if nodeTag ∉ keys(nodeToElementMap)
                    nodeToElementMap[nodeTag] = Vector{Tuple{Tuple{Int64, Int64}, Int64}}(undef, 0)
                end
                Base.push!(nodeToElementMap[nodeTag], (attribute, elementNo))
            end
        end
    end
    return nodeToElementMap
end

function getNodeToElementMap!(meshExtra::MeshExtra, mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1})
    meshExtra.nodeToElementMap = getNodeToElementMap(mesh, attributeArray)
end

function getElementDimMap()
    elementDimMap = Dict{DataType, UnitRange{Int64}}()
    elementDimMap[LineElement] = 1:2
    elementDimMap[TriElement] = 1:3
    elementDimMap[QuadElement] = 1:4
    elementDimMap[TetElement] = 1:4
    elementDimMap[HexElement] = 1:8
    return elementDimMap
end

function getAllFaces(mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1})
    allFaces = Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}()
    attribDim = attributeArray[1][1]
    elementDimMap = getElementDimMap()
    for attribute ∈ attributeArray
        if attribDim == attribute[1]
            elementNo = 1
            for element ∈ mesh.Elements[attribute...]
                for face ∈ collect(combinations(element.nodeTags[elementDimMap[typeof(element)]], attribDim))
                    sortedFace = sort(face)
                    if sortedFace ∉ keys(allFaces)
                        allFaces[sortedFace] = Vector{Tuple{Tuple{Int64, Int64}, Int64}}(undef, 0)
                    end
                    push!(allFaces[sortedFace], (attribute, elementNo))
                end
                elementNo += 1 
            end
        else
            error("Only one dimension type in attribute[1] is allowed.")
        end
    end
    return allFaces
end

"""Returns a list of all the faces that are not shared by any other face, hence the boundaries.
Does not include boundaries of two different attributes or materials."""
function getBoundaryFaces(allFaces::Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}})
    return findall(x->length(x)==1, allFaces)
end

function getBoundaryFaces(meshExtra::MeshExtra)
    return getBoundaryFaces(meshExtra.allFaces)
end

"""Returns a map(Dict) of Boundary face to the  (atrribute1, elementNumberWithAttribute1, atrribute2, elementNumberWithAttribute2)"""
function getAllMaterialBoundaryFaces(meshExtra::MeshExtra)
    commonInternalFaces = findall(x->length(x)==2, meshExtra.allFaces)
    materialBoundaryFaces =  Dict{Array{Int64, 1},Tuple{Tuple{Int64, Int64}, Int64, Tuple{Int64, Int64}, Int64}}()
    for commonInternalFace ∈ commonInternalFaces
        #if the two element attributes are not the same they are material boundaries
        attrib1 = meshExtra.allFaces[commonInternalFace][1][1]
        attrib2 = meshExtra.allFaces[commonInternalFace][2][1]
        #println("attrib1 = $attrib1, attrib2 = $attrib2")
        if attrib1 != attrib2
            element1 = meshExtra.allFaces[commonInternalFace][1][2]
            element2 = meshExtra.allFaces[commonInternalFace][2][2]
            materialBoundaryFaces[commonInternalFace] = (attrib1, element1, attrib2, element2)
        end
    end
    return materialBoundaryFaces
end

function getAllBoundaryNodes(boundaryFaces::Array{Array{Int64, 1}, 1})
    #boundaryFaces = findall(x->length(x)==1, allFaces)
    return unique(sort(vcat(boundaryFaces...)))
end

function getAllInternalNodes(boundaryNodes::Array{Int64, 1}, mesh::Mesh)
    return sort(collect(setdiff(keys(mesh.Nodes), boundaryNodes)))
end

"""Returns a map(Dict) of Boundary face to the  (atrribute, elementNumberWithAttribute)"""
function getElNodeTagsToElementMap(mesh::Mesh, attribArray::Vector{Tuple{Int64, Int64}})

    elNodeTagsToElMap = Dict{Vector{Int64}, Tuple{Tuple{Int64, Int64}, Int64}}()
    for attrib ∈ attribArray
        elements = mesh.Elements[attrib]
        for element ∈ elements
            nodeTags = sort(element.nodeTags)
            elNodeTagsToElMap[nodeTags] = (attrib, nodeTags)
        end
    end
    return elNodeTagsToElMap
end

