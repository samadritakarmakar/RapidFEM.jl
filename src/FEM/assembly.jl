
function assembleMatrix(attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}, mesh::Mesh, localMatrixFunc::Function, problemDim::Int64)
    K_COO::SparseMatrixCOO = SparseMatrixCOO()
    #K_local::Array{Float64,2} = Array{Float64,2}(undef, 0, 0)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArray::Array{Float64,2} = getCoordArray(mesh, element)
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, lagrange)
        K_local::Array{Float64,2} = localMatrixFunc(problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        FEMSparse.assemble_local_matrix!(K_COO, vNodes, vNodes, K_local)
    end
    return SparseArrays.SparseMatrixCSC(K_COO)
end

function assembleVector(problemfunction::Function, attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}, mesh::Mesh, localVectorFunc::Function, problemDim::Int64)::Vector
    f::Vector = zeros(mesh.noOfNodes*problemDim)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArray::Array{Float64,2} = getCoordArray(mesh, element)
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, lagrange)
        f_local::Array{Float64,1} = localVectorFunc(problemfunction, problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        f[vNodes] = f_local
    end
    return f
end
