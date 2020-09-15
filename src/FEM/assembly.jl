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
    numOfThreads::Int64 = Threads.nthreads()    #Total number of threads running
    #An array of SparseMatrixCOO is used to addup matrices of in each thread
    K_COO::Array{SparseMatrixCOO, 1} = Array{SparseMatrixCOO, 1}(undef, numOfThreads)
    Threads.@threads for thread ∈ 1:numOfThreads
        K_COO[thread] = SparseMatrixCOO() #Intiailization of SparseMatrixCOO
    end
    #Intiailization of FeSpaceThreaded
    FeSpaceThreaded::Array{Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, 1} =
    Array{Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, 1}(undef, numOfThreads)
    FeSpaceThreaded[1] = FeSpace #First thread can be the same as the original FeSpace
    Threads.@threads for thread ∈ 2:numOfThreads
        FeSpaceThreaded[thread] = deepcopy(FeSpace) #Others need to be deepcpoies
    end
    RangeDict = createDimRange() # Creates a dim range. Useful when using for 1D or 2D problems
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    #######This makes sure the full size of K_COO is same as the size of the filled K matrix
    vNodes2::Array{Int64,1} = [mesh.noOfNodes*problemDim]
    K_local2::Array{Float64,2} = zeros(1,1)
    ####
    Threads.@threads for thread ∈ 1:numOfThreads #Setting the matrices size for each thread
        FEMSparse.assemble_local_matrix!(K_COO[thread], vNodes2, vNodes2, K_local2)
    end
    ### Multi-Thread Assembly starts here
    Threads.@threads for elementNo ∈ 1:length(mesh.Elements[attribute])
        currentThread::Int64 = Threads.threadid()
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, lagrange)
        K_local::Array{Float64,2} = localMatrixFunc(parameterFunction, problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        FEMSparse.assemble_local_matrix!(K_COO[currentThread], vNodes, vNodes, K_local)
    end
    ### Multi-Thread Assembly ends here
    for thread ∈ 2:numOfThreads ### Adding all the thread matrices to 1st K_COO
        append!(K_COO[1].I, K_COO[thread].I)
        append!(K_COO[1].J, K_COO[thread].J)
        append!(K_COO[1].V, K_COO[thread].V)
        merge!(FeSpace, FeSpaceThreaded[thread]) ##Merge all threaded FeSpace together
    end
    return SparseArrays.SparseMatrixCSC(K_COO[1])
end

function assembleVector(problemfunction::T, attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, mesh::Mesh, localVectorFunc::Function, problemDim::Int64, activeDimensions::Array{Int64,1}=[1, 1, 1])::Vector where T
    numOfThreads::Int64 = Threads.nthreads()    #Total number of threads running
    f::Array{Vector,1} = Array{Vector,1}(undef, numOfThreads)
    Threads.@threads for thread ∈ 1:numOfThreads
        f[thread] = zeros(mesh.noOfNodes*problemDim)
    end
    #Intiailization of FeSpaceThreaded
    FeSpaceThreaded::Array{Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, 1} =
    Array{Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, 1}(undef, numOfThreads)
    FeSpaceThreaded[1] = FeSpace #First thread can be the same as the original FeSpace
    Threads.@threads for thread ∈ 2:numOfThreads
        FeSpaceThreaded[thread] = deepcopy(FeSpace) #Others need to be deepcpoies
    end
    RangeDict = createDimRange()
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    Threads.@threads for elementNo ∈ 1:length(mesh.Elements[attribute])
        currentThread::Int64 = Threads.threadid()
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, lagrange)
        f_local::Array{Float64,1} = localVectorFunc(problemfunction, problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        f[currentThread][vNodes] += f_local
    end
    for thread ∈ 2:numOfThreads ### Adding all the thread vectors to 1st f
        f[1] += f[thread]
        merge!(FeSpace, FeSpaceThreaded[thread]) ##Merge all threaded FeSpace together
    end
    return f[1]
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
