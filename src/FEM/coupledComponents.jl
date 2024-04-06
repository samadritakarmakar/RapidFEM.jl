"""CoupledComponents is used in a global sense and is useful if coupled problems are being used. 
The problem may also be a mixed finite element problem. The struct stores the global end position of solution.
c[i], where is a CoupledComponents variable and i the i_th solution variable, outputs a UnitRange m:n.
Here m is the starting global index and n is the end global index of the given solution variable[i] in 
a global FEM vector
Same can be done for a global FEM matrix c[i,j] which outputs (m_i:n_i,m_j:n_j) where i and j are the 
i_th and j_th respectively solution variables.
"""

struct CoupledComponents
    dofs::Array{Int64, 1}
end

function CoupledComponents(mesh::Mesh, problemDims::Array{Int64, 1})

    dofs = zeros(Int64, length(problemDims))
    lastDof = 0
    for i ∈ 1:length(dofs)
        dofs[i] = mesh.noOfNodes*problemDims[i]+lastDof
        lastDof = dofs[i]
    end
    return CoupledComponents(dofs)
end

function CoupledComponents(meshTuple::NTuple{T1, Mesh}, problemDims::NTuple{T2, Int64}, meshNos::NTuple{T3, Int64}) where {T1, T2, T3}
    @assert length(problemDims) == length(meshNos) "Length of problemDims must be the same as meshNos."
    @assert maximum(meshNos) <= length(meshTuple) "Maximum of meshNos cannot be greater than the total number of meshes."

    dofs = zeros(Int64, length(problemDims))
    lastDof = 0
    for i ∈ 1:length(dofs)
        dofs[i] = meshTuple[meshNos[i]].noOfNodes*problemDims[i]+lastDof
        lastDof = dofs[i]
    end
    return CoupledComponents(dofs)
end


function Base.getindex(c::CoupledComponents, i::Integer)
    lastIndex = i == 1 ? 1 : c.dofs[i-1]+1
    @assert i <=length(c.dofs) "Attempt to access $(length(c.dofs))-length at index $i."
    return lastIndex:c.dofs[i]
end

function Base.getindex(c::CoupledComponents, i::Integer, j::Integer)
    return c[i], c[j]
end

"""Creates a global sparse zero matrix from the CoupledComponents struct variable"""
function SparseArrays.spzeros(c::CoupledComponents)
    return spzeros(c.dofs[end], c.dofs[end])
end

"""Creates a global dense zero vector from the CoupledComponents struct variable"""
function Base.zeros(c::CoupledComponents)
    return zeros(c.dofs[end])
end
