
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
    #=NodeList::Array{Int64} = []
    for element ∈ mesh.Elements[attribute...]
        nodes::Array{Int64} = getNodes(element)
        for node ∈ nodes
            if node ∉ (NodeList)
                push!(NodeList, node)
            end
        end
    end=#
    NodeList::Array{Int64} = []
    for element ∈ mesh.Elements[attribute...]
        nodes::Array{Int64} = getNodes(element)
        for node ∈ nodes
            if length(searchsorted(NodeList, node))==0
                push!(NodeList, node)
                sort!(NodeList)
            end
        end
    end
    #unique!(NodeList) #Makes sure of only one copy of each node in NodeList
    return NodeList
end


"""Applies Dirichlet Boundary condition on the given Stiffness matrix 'A' and
vector 'b' as per the given DirichletFunction which is depedent on the position, x

    applyDirichletBC!(A, b, DirichletFunction, attribute, mesh, problemDim)
"""
function applyDirichletBC!(b::Vector, A::SparseMatrixCSC,
    DirichletFunction::Function, attribute::Tuple{Int64, Int64}, mesh::Mesh,
    problemDim::Int64)

    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim)
    for nodeNo ∈ 1:length(nodes)
        coordArray::Array{Float64} = mesh.Nodes[nodes[nodeNo]]
        b[vNodes[nodeNo:nodeNo+problemDim-1]] = DirichletFunction(coordArray)
    end

    #=A[:, vNodes] .= 0.0
    A[vNodes, :] .= 0.0
    Threads.@threads for vNode ∈ vNodes
        A[vNode, vNode] = 1.0
    end
    return nothing
    =#
    P::SparseMatrixCSC, Pzeros::SparseMatrixCSC =  getPermutionMatrix(vNodes, mesh, problemDim)
    A = P'*A
    A *= P
    A += Pzeros
    return A

end
