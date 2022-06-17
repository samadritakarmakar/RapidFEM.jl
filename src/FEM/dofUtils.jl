#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
"""This function returns the node tags for the given element.

    nodes::Array{Int64} = getNodes(element)
"""
function getNodes(element::AbstractElement)::Array{Int64}
    return element.nodeTags
end

function getAppliedDofLength(problemDim::Int64,
    appliedDof::Array{Int64, 1} = ones(Int, problemDim))
    lengthAppldDof::Int64 = 0
    #println("appliedDof = ", appliedDof, )
    for i ∈ 1:problemDim
        if appliedDof[i] == 1
            lengthAppldDof +=1
        end
    end
    #println("lengthAppldDof = ", lengthAppldDof)
    return lengthAppldDof
end

"""This function returns the vector nodes or the vector degrees of freedom for the
given nodes. This data is important in determining the postion of the data in the
global stiffness matrix and vector.

    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim)
    vNodes::Array{Int64} = getVectorNodes(element, problemDim)
"""
function getVectorNodes(nodes::Array{Int64}, problemDim::Int64,
    appliedDof::Array{Int64, 1} = ones(Int, problemDim))::Array{Int64}

    lengthAppldDof::Int64 = getAppliedDofLength(problemDim, appliedDof)

    vectorNodes::Array{Int64} = Array{Int64}(undef, length(nodes)*lengthAppldDof)
    for i ∈ 1:length(nodes)
        for j ∈ 1:problemDim
            if (appliedDof[j] == 1)
                #sum(appliedDof[1:j])] is a trick to give the current dof number 
                vectorNodes[lengthAppldDof*(i-1)+sum(appliedDof[1:j])] = problemDim*(nodes[i]-1)+j
            end
        end
    end
    return vectorNodes
end
"""This function returns the vector nodes or the vector degrees of freedom for the
given element. This data is important in determining the postion of the data in the
global stiffness matrix and vector.

    vNodes::Array{Int64} = getVectorNodes(nodes, problemDim)
    vNodes::Array{Int64} = getVectorNodes(element, problemDim)
"""
function getVectorNodes(element::AbstractElement, problemDim::Int64,
    appliedDof::Array{Int64, 1} = ones(Int, problemDim))::Array{Int64}

    nodes::Array{Int64} = getNodes(element)
    return getVectorNodes(nodes, problemDim, appliedDof)
end


"""Extracts the solution available at a particular Element for a certain problem dimension.

    solAtNodes::Array{Float64,1} = getSolAtElement(sol, element, problemDim)
"""
function getSolAtElement(sol::AbstractArray{Float64,1}, element::AbstractElement, problemDim::Int)::Array{Float64,1}
    vectorNodes::Array{Int64,1} = getVectorNodes(element, problemDim)
    #sort!(vectorNodes)
    solAtNodes::Array{Float64,1} = Array{Float64,1}(undef, length(vectorNodes))
    i::Int64 = 1
    for node ∈ vectorNodes
        solAtNodes[i] = sol[node]
        i +=1
    end
    return solAtNodes
end

"""Extracts the solution available at a particular Element for a certain problem dimension. 
activeDimensions can be used if the solution dimension is not the same as used processing dimension.
For example, the displacement solution is in 2D but the stesses usually need displacement in 3D.

    solAtNodes::Array{Float64,1} = getSolAtElement(sol, element, problemDim, activeDimensions)
"""
function getSolAtElement(sol::AbstractArray{Float64,1}, element::AbstractElement, problemDim::Int, 
    activeDimensions::Vector{Int64})::Array{Float64,1}

    vectorNodes::Array{Int64,1} = getVectorNodes(element, problemDim)
    #sort!(vectorNodes)
    solAtNodes::Array{Float64,1} = zeros(length(activeDimensions)*length(element.nodeTags))
    j::Int64 = 1
    dimNo = 1 #keeps track of current dim in activeDimensions
    for i ∈ 1:length(solAtNodes)
        if activeDimensions[dimNo] == 0
            solAtNodes[i] = 0.0
        else
            node = vectorNodes[j]
            solAtNodes[i] = sol[node]
            j += 1
        end
        #if dimNo has reached the end then set it to zero
        if dimNo == sum(activeDimensions)
            dimNo = 1
        else
            dimNo += 1
        end
    end
    return solAtNodes
end

"""Extracts views of global solutions in coupled problems from a single mesh.

In the below example u has dimension 3, p has 1 and t has 1.

        u, p, t = getCoupledGlobalSols(globalSol, mesh, [3, 1, 1])
"""
function getCoupledGlobalSols(globalSol::AbstractArray{Float64, 1}, mesh::Mesh, problemDims::Array{Int64, 1})
    sols = Array{AbstractArray{Float64, 1}, 1}(undef, length(problemDims))
    lastSolIndex = 0
    for i ∈ 1:length(sols)
        startIndex = 1+lastSolIndex
        endIndex = problemDims[i]*mesh.noOfNodes+lastSolIndex
        lastSolIndex = endIndex
        sols[i] = @view globalSol[startIndex:endIndex]
    end
    return sols
end

function getCoupledGlobalSols(globalSol::AbstractArray{Float64, 1}, coupling::CoupledComponents)
    sols = Array{AbstractArray{Float64, 1}, 1}(undef, length(coupling.dofs))
    for i ∈ 1:length(sols)
        sols[i] = @view globalSol[coupling[i]]
    end
    return sols
end

function getCoupledGlobalSols(globalSol::AbstractArray{Float64, 1}, meshTuple::Tuple, problemDims::Array{Int64, 1}, meshNos::Array{Int64, 1})
    @assert length(problemDims) == length(meshNos) "Length of problemDims must be the same as meshNos."
    @assert maximum(meshNos) <= length(meshTuple) "Maximum of meshNos cannot be greater than the total number of meshes."
    sols = Array{AbstractArray{Float64, 1}, 1}(undef, length(problemDims))
    lastSolIndex = 0
    for i ∈ 1:length(sols)
        startIndex = 1+lastSolIndex
        endIndex = problemDims[i]*meshTuple[meshNos[i]].noOfNodes+lastSolIndex
        lastSolIndex = endIndex
        sols[i] = @view globalSol[startIndex:endIndex]
    end
    return sols
end