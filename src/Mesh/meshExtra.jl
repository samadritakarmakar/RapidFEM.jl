using Combinatorics

mutable struct MeshExtra
    nodeToElementMap::Dict{Int64, Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}}

    
    #MeshExtra(nodeToElementMap::Dict{Int64, Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}}) = new(nodeToElementMap)
end

MeshExtra() = MeshExtra(Dict{Int64, Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}}())

MeshExtra(mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1}) = MeshExtra(getNodeToElementMap(mesh, attributeArray))

function getNodeToElementMap(mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64}, 1})
    
    nodeToElementMap = Dict{Int64, Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}}()
    for attribute ∈ attributeArray
        elements = mesh.Elements[attribute...]
        noOfElements = size(elements,1)
        for elementNo ∈ 1:noOfElements
            element = elements[elementNo]
            for nodeTag ∈ element.nodeTags
            #for nodeNo ∈ 1:length(element)
                #nodeTag = element.nodeTags[nodeNo]
                if nodeTag ∉ keys(nodeToElementMap)
                    nodeToElementMap[nodeTag] = Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}(undef, 0)
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
    allFaces = Dict{Array{Int64, 1}, Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}}()
    attribDim = attributeArray[1][1]
    elementDimMap = getElementDimMap()
    for attribute ∈ attributeArray
        if attribDim == attribute[1]
            elementNo = 1
            for element ∈ mesh.Elements[attribute...]
                for face ∈ collect(combinations(element.nodeTags[elementDimMap[typeof(element)]], attribDim))
                    if face ∉ keys(allFaces)
                        allFaces[face] = Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}(undef, 0)
                    end
                    push!(allFaces[face], (attribute, elementNo))
                end
                elementNo += 1 
            end
        end
    end
    return allFaces
end

function getAllBoundaryNodes(allFaces::Dict{Array{Int64, 1}, Array{Tuple{Tuple{Int64, Int64}, Int64}, 1}})
    boundaryFaces = findall(x->length(x)==1, allFaces)
    return unique(sort(vcat(boundaryFaces...)))
end

function getAllInternalNodes(boundaryNodes::Array{Int64}, mesh::Mesh)
    return setdiff(keys(mesh.Nodes), boundaryNodes)    
end