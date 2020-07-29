function getNodes(element::AbstractElement)::Array{Int64}
    return element.nodeTags
end

function getVectorNodes(nodes::Array{Int64}, problemDim::Int64)::Array{Int64}
    vectorNodes::Array{Int64} = Array{Int64}(undef, length(nodes)*problemDim)
    for i ∈ 1:length(nodes)
        for j ∈ 1:problemDim
            vectorNodes[problemDim*(i-1)+j] = problemDim*(nodes[i]-1)+j
        end
    end
    return vectorNodes
end

function getVectorNodes(element::AbstractElement, problemDim::Int64)::Array{Int64}
    nodes::Array{Int64} = getNodes(element)
    return getVectorNodes(nodes, problemDim)
end
