include("meshExtra.jl")

"""Replaces and Adds Elements in a Set of element attibutes. If replaceElAttribsElNos is used directly then it is recommended 
that all element attributes be the same.

    replaceAndAddElements!(mesh, replaceElAttribsElNos, newElementNodeTags)

    replaceAndAddElements!(mesh, attrib, elementNos, newElementNodeTags)
"""

function replaceAndAddElements!(mesh::Mesh, 
    replaceElAttribsElNos::Set{Tuple{Tuple{Int64, Int64}, Int64}},
    newElementNodeTags::Vector{Vector{Int64}})

    newElAttribsElNos = Vector{Tuple{Tuple{Int64, Int64}, Int64}}()

    lengthNewElementTags = length(newElementNodeTags)
    newElementTagNo = 1
    deletedElements = 0
    attrib = (0,0)
    for replaceElAttribsElNo ∈ replaceElAttribsElNos
        attrib, changeTagElNo = replaceElAttribsElNo
        if newElementTagNo <= lengthNewElementTags
            label = mesh.Elements[attrib][changeTagElNo].label
            element = createNewElement(attrib, label, newElementNodeTags[newElementTagNo])
            
            mesh.Elements[attrib][changeTagElNo] = element
            
            push!(newElAttribsElNos, (attrib, changeTagElNo))

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
        push!(newElAttribsElNos, length(mesh.Elements[attrib]))   
    end
    mesh.noOfElements += addElements - deletedElements

    return newElAttribsElNos
end

function replaceAndAddElements!(mesh::Mesh, attrib::Tuple{Int64, Int64},
    elementNos::Union{Vector{Int64}, Set{Int64}}, 
    newElementNodeTags::Vector{Vector{Int64}})

    replaceElAttribsElNos = Set{Tuple{Tuple{Int64, Int64}, Int64}}()
    for elementNo ∈ elementNos
        push!(replaceElAttribsElNos, (attrib, elementNo))
    end
    replaceAndAddElements!(mesh, replaceElAttribsElNos, newElementNodeTags)
end

function getBoundaryPolyhedron(elementNodeTagsArray::Vector{Vector{Int64}}, dimension::Int64)
    boundaryPolyhedron = Set{Vector{Int64}}()
    for elementNodeTags ∈ elementNodeTagsArray
        elementType, order = getElementTypeAndOrder(dimension, elementNodeTags)
        if elementType == LineElement ||elementType == PointElement || order > 1
            error("element Type $elementType or order $order not supported to find outer polygons")
        end
        combinations = getFaceCombinations(elementType, order)
        faces = getFacesFromCombinations(elementNodeTags, combinations)
        for face ∈ faces
            if face ∉ boundaryPolyhedron
                push!(boundaryPolyhedron, face)
            else
                pop!(boundaryPolyhedron, face)
            end
        end
    end
    return boundaryPolyhedron
end




