function getPermutionMatrix(vNodes::Array{Int64}, mesh::Mesh, problemDim::Int64)
    noOfVectorNodes::Int64 = mesh.noOfNodes*problemDim
    Pzeros::SparseMatrixCSC = spzeros(noOfVectorNodes, noOfVectorNodes)
    Pzeros[vNodes, vNodes] = sparse(I, length(vNodes), length(vNodes))
    P::SparseMatrixCSC = sparse(I, noOfVectorNodes, noOfVectorNodes)
    P = P - Pzeros
    return P, Pzeros
end

"""Generates an array of Unique nodes for the given attribute in the mesh data."""
function getUniqueNodes(attribute::Tuple{Int64, Int64}, mesh::Mesh)::Array{Int64}
    NodeList::Array{Int64} = []
    for element ∈ mesh.Elements[attribute...]
        nodes::Array{Int64} = getNodes(element)
        for node ∈ nodes
            if node ∉ (NodeList)
                push!(NodeList, node)
            end
        end
    end
    return NodeList
end

function applyDirichletBC!(A::SparseMatrixCSC, b::Vector, DirichletFunction::Function, attribute::Tuple{Int64, Int64}, mesh::Mesh, problemDim::Int64)
    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim)
    for nodeNo ∈ 1:length(nodes)
        coordArray::Array{Float64} = mesh.Nodes[nodes[nodeNo]]
        b[vNodes[nodeNo:nodeNo+problemDim-1]] = DirichletFunction(coordArray)
    end
    P::SparseMatrixCSC, Pzeros::SparseMatrixCSC =  getPermutionMatrix(vNodes, mesh, problemDim)
    A[:] = P'*A*P
    A[:] = A + Pzeros
end
