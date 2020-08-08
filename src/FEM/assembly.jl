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

function assembleMatrix(parameterFunction::T, attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, mesh::Mesh, localMatrixFunc::Function, problemDim::Int64, activeDimensions::Array{Int64,1}=[1, 1, 1])::SparseMatrixCSC where T
    K_COO::SparseMatrixCOO = SparseMatrixCOO()
    RangeDict = createDimRange()
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    #K_local::Array{Float64,2} = Array{Float64,2}(undef, 0, 0)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, lagrange)
        K_local::Array{Float64,2} = localMatrixFunc(parameterFunction, problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        FEMSparse.assemble_local_matrix!(K_COO, vNodes, vNodes, K_local)
    end
    vNodes2::Array{Int64,1} = [mesh.noOfNodes*problemDim]
    K_local2::Array{Float64,2} = zeros(1,1)
    FEMSparse.assemble_local_matrix!(K_COO, vNodes2, vNodes2, K_local2)
    return SparseArrays.SparseMatrixCSC(K_COO)
end

function assembleVector(problemfunction::T, attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, mesh::Mesh, localVectorFunc::Function, problemDim::Int64, activeDimensions::Array{Int64,1}=[1, 1, 1])::Vector where T
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
        f[vNodes] += f_local
    end
    return f
end


function assembleScalar(problemfunction::T, attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, mesh::Mesh, localVectorFunc::Function, problemDim::Int64, activeDimensions::Array{Int64,1}=[1, 1, 1])::Vector where T
    noOfElements::Int64 = getNoOfElements(mesh, attribute)
    f::Vector = zeros(noOfElements*problemDim)
    RangeDict = createDimRange()
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, lagrange)
        f_local::Array{Float64,1} = localVectorFunc(problemfunction, problemDim, element, shapeFunction, coordArray)
        f[problemDim*(elementNo-1)+1:problemDim*elementNo] += f_local
    end
    return f
end
