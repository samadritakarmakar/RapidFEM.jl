struct CoupledComponents
    dofs::Array{Int64, 1}
    function CoupledComponents(mesh::Mesh, problemDims::Array{Int64, 1})

        dofs = zeros(Int64, length(problemDims))
        lastDof = 0
        for i ∈ 1:length(dofs)
            dofs[i] = mesh.noOfNodes*problemDims[i]+lastDof
            lastDof = dofs[i]
        end
        return new(dofs)
    end

    function CoupledComponents(meshTuple::Tuple, problemDims::Array{Int64, 1}, meshNos::Array{Int64, 1})
        @assert length(problemDims) == length(meshNos) "Length of problemDims must be the same as meshNos."
        @assert maximum(meshNos) <= length(meshTuple) "Maximum of meshNos cannot be greater than the total number of meshes."

        dofs = zeros(Int64, length(problemDims))
        lastDof = 0
        for i ∈ 1:length(dofs)
            dofs[i] = meshTuple[meshNos[i]].noOfNodes*problemDims[i]+lastDof
            lastDof = dofs[i]
        end
        return new(dofs)
    end
end

function Base.getindex(c::CoupledComponents, i::Integer)
    lastIndex = i == 1 ? 1 : c.dofs[i-1]+1
    @assert i <=length(c.dofs) "Attempt to access $(length(c.dofs))-length at index $i."
    return lastIndex:c.dofs[i]
end

function Base.getindex(c::CoupledComponents, i::Integer, j::Integer)
    return c[i], c[j]
end