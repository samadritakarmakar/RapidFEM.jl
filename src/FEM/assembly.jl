function createDimRange()::Dict{Array{Int64,1}, StepRange{Int64,Int64}}
    RangeDict::Dict{Array{Int64,1}, StepRange{Int64,Int64}} = Dict()
    if length(keys(RangeDict)) < 7
        RangeDict[[1, 1, 1]] = 1:1:3
        RangeDict[[0, 1, 1]] = 2:1:3
        RangeDict[[1, 0, 1]] = 1:2:3
        RangeDict[[1, 1, 0]] = 1:1:2
        RangeDict[[0, 0, 1]] = 3:1:3
        RangeDict[[0, 1, 0]] = 2:1:2
        RangeDict[[1, 0, 0]] = 1:1:1
    end
    return RangeDict
end

function getRange(RangeDict::Dict{Array{Int64,1}, StepRange{Int64,Int64}}, activeDimensions::Array{Int64,1})::StepRange{Int64,Int64}
    return RangeDict[activeDimensions]
end

function assembleMatrix(attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}, mesh::Mesh, localMatrixFunc::Function, problemDim::Int64, activeDimensions::Array{Int64,1}=[1, 1, 1])::SparseMatrixCSC
    K_COO::SparseMatrixCOO = SparseMatrixCOO()
    RangeDict = createDimRange()
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    #K_local::Array{Float64,2} = Array{Float64,2}(undef, 0, 0)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, lagrange)
        K_local::Array{Float64,2} = localMatrixFunc(problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        FEMSparse.assemble_local_matrix!(K_COO, vNodes, vNodes, K_local)
    end
    return SparseArrays.SparseMatrixCSC(K_COO)
end

function assembleVector(problemfunction::Function, attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}, mesh::Mesh, localVectorFunc::Function, problemDim::Int64, activeDimensions::Array{Int64,1}=[1, 1, 1])::Vector
    f::Vector = zeros(mesh.noOfNodes*problemDim)
    RangeDict = createDimRange()
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, lagrange)
        f_local::Array{Float64,1} = localVectorFunc(problemfunction, problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        f[vNodes] = f_local
    end
    return f
end
