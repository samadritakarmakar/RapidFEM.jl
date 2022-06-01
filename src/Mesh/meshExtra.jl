using Combinatorics

struct MeshExtra
    nodeToElementMap::Dict{Int64, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}
    allFaces::Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}
    #boundaryFaces::Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}
    boundaryNodes::Array{Int64, 1}
    internalNodes::Array{Int64, 1}
    #elementNodeTagsToElementMap::Dict{Vector{Int64}, Tuple{Tuple{Int64, Int64}, Int64}}
    #MeshExtra(nodeToElementMap::Dict{Int64, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}) = new(nodeToElementMap)
end

function MeshExtra()
    MeshExtra(Dict{Int64, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}(), 
    Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}(),
    zeros(Int64, 0), zeros(Int64, 0))#=, 
    Dict{Vector{Int64}, Tuple{Tuple{Int64, Int64}, Int64}}())=#
end

function MeshExtra(mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1})
    nodeToElementMap = getNodeToElementMap(mesh, attributeArray)
    allFaces = getAllFaces(mesh, attributeArray)
    boundaryFaces = getBoundaryFaces(allFaces)
    boundaryNodes = getAllBoundaryNodes(boundaryFaces)
    internalNodes = getAllInternalNodes(boundaryNodes, mesh)
    #elNodeTagsToElMap = getElNodeTagsToElementMap(mesh::Mesh, attributeArray)
    MeshExtra(nodeToElementMap, allFaces, boundaryNodes, internalNodes)#, elNodeTagsToElMap)
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

function getFacesFromCombinations(nodes::Vector{Int64}, combinations::Vector{Vector{Int64}})
    faces = zero.(combinations)
    combinationNo = 1
    for combination ∈ combinations
        faces[combinationNo] .= sort(nodes[combination])
        combinationNo += 1
    end
    return faces
end

function getFaceCombinations(elementType::DataType, order::Int64)
    if elementType == LineElement
        combinations = [[1], [2]]
        if order == 2
            combinations = [[1], [3], [2]]
        elseif order == 3
            combinations = [[1], [3], [4], [2]]
        elseif order != 1
            error("Order $(order) not supported.")
        end
    if elementType == TriElement
        combinations = [[1, 2], [1, 3], [2, 3]]
        if order == 2
            combinations = [[1, 2, 4], [1, 3, 6], [2, 3, 5]]
        elseif order == 3
            combinations = [[1, 2, 4, 5], [1, 3, 8, 9], [2, 3, 6, 7]]
        elseif order != 1
            error("Order $(order) not supported.")
        end
    elseif elementType == QuadElement
        combinations = [[1, 2], [2, 3], [3, 4], [1, 4]]
        if order == 2
            combinations = [[1, 2, 5], [2, 3, 6], [3, 4, 7], [1, 4, 8]]
        elseif order == 3
            combinations = [[1, 2, 5, 6], [2, 3, 7, 8], [3, 4, 9, 10], [1, 4, 11, 12]]
        elseif order != 1
            error("Order $(order) not supported.")
        end
    elseif elementType == TetElement
        combinations = [[1, 2, 4], [1, 2, 3],[1, 3, 4], [2, 3, 4]]
        if order == 2
            combinations = [[1, 2, 4, 5, 8, 10], [1, 2, 3, 5, 6, 7], [1, 3, 4, 7, 8, 9], [2, 3, 4, 6, 9, 10]]
        elseif order == 3
            combinations = [[1, 2, 4, 5, 6, 11, 12, 15, 16, 18], 
            [1, 3, 4, 9, 10, 11, 12, 13, 14, 19], 
            [2, 3, 4, 7, 8, 13, 14, 15, 16, 20],
            [1, 2, 3, 5, 6, 7, 8, 9, 10, 17]]
        elseif order != 1
            error("Order $(order) not supported.")
        end
    elseif elementType == HexElement
        combinations = [[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 5, 6], [3, 4, 7, 8], [1, 4, 5, 8], [2, 3, 6, 7]]
        if order == 2
            combinations = [[1, 2, 3, 4, 9, 10, 12, 14, 21], [2, 3, 6, 7, 12, 13, 15, 19, 24], 
            [1, 4, 5, 8, 10, 11, 16, 18, 23], [1, 2, 5, 6, 9, 11, 13, 17, 22],
            [5, 6, 7, 8, 17, 18, 19, 20, 26], [3, 4, 7, 8, 14, 15, 16, 20, 25]]
        elseif order == 3
            combinations = [[1, 4, 5, 8, 11, 12, 13, 14, 23, 24, 27, 28, 41, 42, 43, 44], 
            [2, 3, 6, 7, 15, 16, 17, 18, 21, 22, 29, 30, 45, 46, 47, 48],
            [1, 2, 3, 4, 9, 10, 11, 12, 15, 16, 19, 20, 33, 34, 35, 36],
            [5, 6, 7, 8, 25, 26, 27, 28, 29, 30, 31, 32, 53, 54, 55, 56],
            [1, 2, 5, 6, 9, 10, 13, 14, 17, 18, 25, 26, 37, 38, 39, 40],
            [3, 4, 7, 8, 19, 20, 21, 22, 23, 24, 31, 32, 49, 50, 51, 52]]
        elseif order != 1
            error("Order $(order) not supported.")
        end
    end
    return combinations
end

function getElementFaces(element::AbstractElement)
    combinations = getFaceCombinations(typeof(element), element.order)
    return getFacesFromCombinations(element.nodeTags, combinations)
end

function getElementFaces(element::HexElement)
    combinations = [[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 5, 6], [3, 4, 7, 8], [1, 4, 5, 8], [2, 3, 6, 7]]
    if element.order == 2
        combinations = [[1, 2, 3, 4, 9, 10, 12, 14, 21], [2, 3, 6, 7, 12, 13, 15, 19, 24], 
        [1, 4, 5, 8, 10, 11, 16, 18, 23], [1, 2, 5, 6, 9, 11, 13, 17, 22],
        [5, 6, 7, 8, 17, 18, 19, 20, 26], [3, 4, 7, 8, 14, 15, 16, 20, 25]]
    elseif element.order == 3
        combinations = [[1, 4, 5, 8, 11, 12, 13, 14, 23, 24, 27, 28, 41, 42, 43, 44], 
        [2, 3, 6, 7, 15, 16, 17, 18, 21, 22, 29, 30, 45, 46, 47, 48],
        [1, 2, 3, 4, 9, 10, 11, 12, 15, 16, 19, 20, 33, 34, 35, 36],
        [5, 6, 7, 8, 25, 26, 27, 28, 29, 30, 31, 32, 53, 54, 55, 56],
        [1, 2, 5, 6, 9, 10, 13, 14, 17, 18, 25, 26, 37, 38, 39, 40],
        [3, 4, 7, 8, 19, 20, 21, 22, 23, 24, 31, 32, 49, 50, 51, 52]]
    elseif element.order != 1
        error("Order $(element.order) not supported.")
    end
    return getFacesFromCombinations(element.nodeTags, combinations)
end



function getAllFaces(mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1})
    allFaces = Dict{Array{Int64, 1}, Vector{Tuple{Tuple{Int64, Int64}, Int64}}}()
    attribDim = attributeArray[1][1]
    #elementDimMap = getElementDimMap()
    for attribute ∈ attributeArray
        if attribDim == attribute[1]
            elementNo = 1
            for element ∈ mesh.Elements[attribute...]
                #for face ∈ collect(combinations(element.nodeTags[elementDimMap[typeof(element)]], attribDim))
                faces = getElementFaces(element)
                for face ∈ faces
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

