include("mesh.jl")
include("meshExtra.jl")

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

"""Replaces and Adds Elements in a Set of element attibutes. If replaceElAttribsElNos is used directly then it is recommended 
that all element attributes be the same.

    replaceAndAddElements!(mesh, replaceElAttribsElNos, newElementNodeTags)

    replaceAndAddElements!(mesh, attrib, elementNos, newElementNodeTags)
"""

function replaceAndAddElements!(mesh::Mesh, 
    replaceElAttribsElNos::Set{Tuple{Tuple{Int64, Int64}, Int64}},
    newElementNodeTags::Vector{Vector{Int64}})

    newElAttribsElNos = Vector{Tuple{Tuple{Int64, Int64}, Int64}}()
    if newElementNodeTags isa Set
        newElementNodeTags = collect(newElementNodeTags)
    end

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

function updateElNodeTagsToElementMap!(elNodeTagsToElementMap::Dict{Vector{Int64}, Tuple{Tuple{Int64, Int64}, Int64}},
    replaceElAttribsElNos::Set{Tuple{Tuple{Int64, Int64}, Int64}}, newElAttribsElNos::Vector{Tuple{Tuple{Int64, Int64}, Int64}},
    newElementTagSet::Vector{Vector{Int64}}, mesh::Mesh)


end

function addNewElNodeTagsToElementMap!(elNodeTagsToElementMap::Dict{Vector{Int64}, Tuple{Tuple{Int64, Int64}, Int64}},
    newElementTagSet::Vector{Vector{Int64}}, attrib::Tuple{Int64, Int64})
end


