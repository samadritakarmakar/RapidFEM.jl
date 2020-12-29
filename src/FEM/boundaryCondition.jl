#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

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
    problemDim::Int64,
    appliedDof::Array{Int64, 1} = ones(Int, problemDim), varArgs...)

    lengthAppldDof::Int64 = getAppliedDofLength(problemDim, appliedDof)
    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim, appliedDof)
    for nodeNo ∈ 1:length(nodes)
        coordArray::Array{Float64} = mesh.Nodes[nodes[nodeNo]]
        diricFunc::Array{Float64,1} = DirichletFunction(coordArray, varArgs...)
        for j ∈ 1:problemDim
            if appliedDof[j] == 1
                b[vNodes[(nodeNo-1)*lengthAppldDof+sum(appliedDof[1:j])]] = diricFunc[j]
                #b[vNodes[(nodeNo-1)*lengthAppldDof+sum(appliedDof[1:j])]]  = 0.0
            end
        end
        #b[vNodes[nodeNo:nodeNo+problemDim-1]] = DirichletFunction(coordArray, varArgs...)
        #b[vNodes[(nodeNo-1)*problemDim+1:nodeNo*problemDim]] = DirichletFunction(coordArray, varArgs...)
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

"""Applies Initial Boundary condition on the given vector 'u' as per the given
    InitialFunction which is depedent on the position, x
"""
function applyInitialBC!(u::Vector, InitialBcFunction::Function,
    attribute::Tuple{Int64, Int64}, mesh::Mesh, problemDim::Int64, varArgs...)

    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim)
    for nodeNo ∈ 1:length(nodes)
        coordArray::Array{Float64} = mesh.Nodes[nodes[nodeNo]]
        u[vNodes[(nodeNo-1)*problemDim+1:nodeNo*problemDim]] .= InitialBcFunction(coordArray, varArgs...)
    end
end

function applyNLDirichletBC_on_J!(J::SparseMatrixCSC,
    attribute::Tuple{Int64, Int64}, mesh::Mesh,
    problemDim::Int64, appliedDof::Array{Int64, 1} = ones(Int, problemDim), varArgs...)

    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim, appliedDof)
    P::SparseMatrixCSC, Pzeros::SparseMatrixCSC =  getPermutionMatrix(vNodes, mesh, problemDim)
    J = P'*J
    J *= P
    J += Pzeros
    return J
end

function applyNLDirichletBC_on_Soln!(Soln::Array{Float64,1},
    DirichletFunction::Function,  attribute::Tuple{Int64, Int64},
    mesh::Mesh, problemDim::Int64,
    appliedDof::Array{Int64, 1} = ones(Int, problemDim), varArgs...)

    lengthAppldDof::Int64 = getAppliedDofLength(problemDim, appliedDof)
    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim, appliedDof)
    for nodeNo ∈ 1:length(nodes)
        coordArray::Array{Float64} = mesh.Nodes[nodes[nodeNo]]
        diricFunc::Array{Float64,1} = DirichletFunction(coordArray, varArgs...)
        for j ∈ 1:problemDim
            if appliedDof[j] == 1
                Soln[vNodes[(nodeNo-1)*lengthAppldDof+sum(appliedDof[1:j])]] = diricFunc[j]
            end
        end
        #Soln[vNodes[(nodeNo-1)*problemDim+1:nodeNo*problemDim]] .= DirichletFunction(coordArray, varArgs...)
    end
end

function applyNLDirichletBC_on_f!(f::Vector,
    attribute::Tuple{Int64, Int64}, mesh::Mesh, problemDim::Int64,
    appliedDof::Array{Int64, 1} = ones(Int, problemDim), varArgs...)

    lengthAppldDof::Int64 = getAppliedDofLength(problemDim, appliedDof)
    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim, appliedDof)
    #println("vNodes = ",vNodes)
    for nodeNo ∈ 1:length(nodes)
        for j ∈ 1:problemDim
            if appliedDof[j] == 1
                f[vNodes[(nodeNo-1)*lengthAppldDof+sum(appliedDof[1:j])]] = 0.0
            end
        end
        #f[vNodes[(nodeNo-1)*problemDim+1:nodeNo*problemDim]] .= 0.0
        #println("Dirichlet vNodes = ", vNodes[(nodeNo-1)*problemDim+1:nodeNo*problemDim])
    end
end

function applyDynamicDirichletBC!(SolutionArray::Array{Array{Float64,1},1},
    b::Vector, A::SparseMatrixCSC, DirichletFunction::Function,
    attribute::Tuple{Int64, Int64}, mesh::Mesh, problemDim::Int64,
    appliedDof::Array{Int64, 1} = ones(Int, problemDim), varArgs...)

    lengthAppldDof::Int64 = getAppliedDofLength(problemDim, appliedDof)
    nodes::Array{Int64} = getUniqueNodes(attribute, mesh)
    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim, appliedDof)
    for nodeNo ∈ 1:length(nodes)
        coordArray::Array{Float64} = mesh.Nodes[nodes[nodeNo]]
        diricFunc::Array{Float64,1} = DirichletFunction(coordArray, varArgs...)
        for j ∈ 1:problemDim
            if appliedDof[j] == 1
                Soln[vNodes[(nodeNo-1)*lengthAppldDof+sum(appliedDof[1:j])]] = diricFunc[j]
                b[vNodes[(nodeNo-1)*lengthAppldDof+sum(appliedDof[1:j])]]  = 0.0
            end
        end
        #SolutionArray[1][vNodes[(nodeNo-1)*problemDim+1:nodeNo*problemDim]] .= DirichletFunction(coordArray, varArgs...)
        #b[vNodes[(nodeNo-1)*problemDim+1:nodeNo*problemDim]] .= 0.0
    end
    P::SparseMatrixCSC, Pzeros::SparseMatrixCSC =  getPermutionMatrix(vNodes, mesh, problemDim)
    A = P'*A
    A *= P
    A += Pzeros
    return A
end
