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
    return getVectorNodes(nodes, problemDim)
end
