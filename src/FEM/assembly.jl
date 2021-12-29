#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
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

"""Responsible for assembling of FEM Matrices. paramerterFunction could
be anything depending on what you would like to send to the local matrix
assembler. The below example is similar to used in: src/Examples/poisson.jl
This function is useful if the local matrix makes changes to the parameters variable.

    K::SparseMatrixCSC = assembleMatrix!(parameters, volAttrib, FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)
"""
function assembleMatrix!(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh,  localMatrixFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1], varArgs...; reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::SparseMatrixCSC where T

    numOfThreads::Int64 = Threads.nthreads()    #Total number of threads running
    #An array of SparseMatrixCOO is used to addup matrices of in each thread
    K_COO::Array{SparseMatrixCOO, 1} = Array{SparseMatrixCOO, 1}(undef, numOfThreads)
    K_localArray::Array{Array{Float64,2},1} = Array{Array{Float64,2},1}(undef, numOfThreads)
    Threads.@threads for thread ∈ 1:numOfThreads
        K_COO[thread] = SparseMatrixCOO() #Intiailization of SparseMatrixCOO
        K_localArray[thread] = zeros(0,0)
    end
    #Intiailization of FeSpaceThreaded
    FeSpaceThreaded::Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1} =
    Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1}(undef, numOfThreads)
    FeSpaceThreaded[1] = FeSpace #First thread can be the same as the original FeSpace
    Threads.@threads for thread ∈ 2:numOfThreads
        FeSpaceThreaded[thread] = deepcopy(FeSpace) #Others need to be deepcopies
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
        #if the type of element changes then reallocate the local matrix else replace with zeros
        if length(K_localArray[currentThread]) != problemDim*element.noOfElementNodes
            K_localArray[currentThread] = zeros(problemDim*element.noOfElementNodes,problemDim*element.noOfElementNodes)
        else
            fill!(K_localArray[currentThread],0.0)
        end
        #K_localArray[currentThread] = zeros(problemDim*element.noOfElementNodes,problemDim*element.noOfElementNodes)
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, reduction = reduction, 
        elementFunction = elementFunction, quadrature = quadrature)
        
        localMatrixFunc(K_localArray[currentThread], parameters, problemDim, element, elementNo, shapeFunction, coordArray, varArgs...)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        FEMSparse.assemble_local_matrix!(K_COO[currentThread], vNodes, vNodes, K_localArray[currentThread])
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

"""Responsible for assembling of FEM Matrices. paramerterFunction could
be anything depending on what you would like to send to the local matrix
assembler. The below example is similar to used in: src/Examples/poisson.jl
This function is useful if the local matrix makes changes to the parameters variable.

    K::SparseMatrixCSC = assembleMatrix!(parameters, volAttrib, FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)
"""
function assembleMatrix(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh,  localMatrixFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1], 
    varArgs...; reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::SparseMatrixCSC where T

    return assembleMatrix!(parameters, attribute,
        FeSpace, mesh, localMatrixFunc, problemDim, activeDimensions, varArgs..., ; 
        reduction = reduction, elementFunction = elementFunction, quadrature = quadrature)
end

"""Responsible for assembling of FEM Vectors. parameters could
be anything depending on what you would like to send to the local matrix
assembler. The below example is similar to used in: src/Examples/poisson.jl
This function is useful if the local matrix makes changes to the parameters variable.

    f::Vector = RapidFEM.assembleVector!(parameters, volAttrib, FeSpace, mesh, RapidFEM.localSource, problemDim, activeDimensions)
"""
function assembleVector!(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh, localVectorFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1],varArgs...; 
    reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::Vector where T

    numOfThreads::Int64 = Threads.nthreads()    #Total number of threads running
    f::Array{Vector,1} = Array{Vector,1}(undef, numOfThreads)
    f_localArray::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, numOfThreads)
    Threads.@threads for thread ∈ 1:numOfThreads
        f[thread] = zeros(mesh.noOfNodes*problemDim)
        f_localArray[thread] = zeros(0)
    end
    #Intiailization of FeSpaceThreaded
    FeSpaceThreaded::Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1} =
    Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1}(undef, numOfThreads)
    FeSpaceThreaded[1] = FeSpace #First thread can be the same as the original FeSpace
    Threads.@threads for thread ∈ 2:numOfThreads
        FeSpaceThreaded[thread] = deepcopy(FeSpace) #Others need to be deepcopies
    end
    RangeDict = createDimRange()
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    Threads.@threads for elementNo ∈ 1:length(mesh.Elements[attribute])
        currentThread::Int64 = Threads.threadid()
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        #if the type of element changes then reallocate the local vector else replace with zeros
        if length(f_localArray[currentThread]) != problemDim*element.noOfElementNodes
            f_localArray[currentThread] = zeros(problemDim*element.noOfElementNodes)
        else
            fill!(f_localArray[currentThread],0.0)
        end
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, reduction = reduction, 
        elementFunction = elementFunction, quadrature = quadrature)

        localVectorFunc(f_localArray[currentThread], parameters, problemDim, element, elementNo, shapeFunction, coordArray, varArgs...)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        f[currentThread][vNodes] += f_localArray[currentThread]
    end
    for thread ∈ 2:numOfThreads ### Adding all the thread vectors to 1st f
        f[1] += f[thread]
        merge!(FeSpace, FeSpaceThreaded[thread]) ##Merge all threaded FeSpace together
    end
    return f[1]
end

"""Responsible for assembling of FEM Vectors. parameters could
be anything depending on what you would like to send to the local matrix
assembler. The below example is similar to used in: src/Examples/poisson.jl

    f::Vector = RapidFEM.assembleVector!(parameters, volAttrib, FeSpace, mesh, RapidFEM.localSource, problemDim, activeDimensions)
"""
function assembleVector(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh, localVectorFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1],varArgs...; 
    reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::Vector where T

    return assembleVector!(parameters, attribute, FeSpace, mesh,
    localVectorFunc, problemDim, activeDimensions, varArgs...; 
    reduction = reduction, elementFunction = elementFunction, quadrature = quadrature)
end


"""Responsible for assembling of FEM Vectors and Matrices with the same domain attribute.
This kind usage is important where the Vectors and Matrices are closely related like in
the case of plasticity
paramerterFunction could be anything depending on what you would like to send to the local matrix
assembler.
This function is useful if the local vector or matrix makes changes to the parameters variable.

    K::SparseMatrixCSC = assembleVectorMatrix!(parameters, volAttrib, FeSpace, mesh, RapidFEM.local_∇v_σ_Vector!,
    RapidFEM.local_∇v_Cᵀ_∇u!, problemDim, activeDimensions)
"""
function assembleVectorMatrix!(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh, localVectorFunc::Function, localMatrixFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1], varArgs...;
    reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss) where T

    numOfThreads::Int64 = Threads.nthreads()    #Total number of threads running
    f::Array{Vector,1} = Array{Vector,1}(undef, numOfThreads)
    f_localArray::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, numOfThreads)
    #An array of SparseMatrixCOO is used to addup matrices of in each thread
    K_COO::Array{SparseMatrixCOO, 1} = Array{SparseMatrixCOO, 1}(undef, numOfThreads)
    K_localArray::Array{Array{Float64,2},1} = Array{Array{Float64,2},1}(undef, numOfThreads)
    Threads.@threads for thread ∈ 1:numOfThreads
        f[thread] = zeros(mesh.noOfNodes*problemDim)
        f_localArray[thread] = zeros(0)
        K_COO[thread] = SparseMatrixCOO() #Intiailization of SparseMatrixCOO
        K_localArray[thread] = zeros(0,0)
    end
    #Intiailization of FeSpaceThreaded
    FeSpaceThreaded::Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1} =
    Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1}(undef, numOfThreads)
    FeSpaceThreaded[1] = FeSpace #First thread can be the same as the original FeSpace
    Threads.@threads for thread ∈ 2:numOfThreads
        FeSpaceThreaded[thread] = deepcopy(FeSpace) #Others need to be deepcopies
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
        #if the type of element changes then reallocate the local vector else replace with zeros
        if length(f_localArray[currentThread]) != problemDim*element.noOfElementNodes
            f_localArray[currentThread] = zeros(problemDim*element.noOfElementNodes)
        else
            fill!(f_localArray[currentThread],0.0)
        end
        #if the type of element changes then reallocate the local matrix else replace with zeros
        if length(K_localArray[currentThread]) != problemDim*element.noOfElementNodes
            K_localArray[currentThread] = zeros(problemDim*element.noOfElementNodes,problemDim*element.noOfElementNodes)
        else
            fill!(K_localArray[currentThread],0.0)
        end
        #K_localArray[currentThread] = zeros(problemDim*element.noOfElementNodes,problemDim*element.noOfElementNodes)
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, reduction = reduction, 
        elementFunction = elementFunction, quadrature = quadrature)

        localVectorFunc(f_localArray[currentThread], parameters, problemDim, element, elementNo, shapeFunction, coordArray, varArgs...)
        localMatrixFunc(K_localArray[currentThread], parameters, problemDim, element, elementNo, shapeFunction, coordArray, varArgs...)
        vNodes::Array{Int64} = getVectorNodes(element, problemDim)
        f[currentThread][vNodes] += f_localArray[currentThread]
        FEMSparse.assemble_local_matrix!(K_COO[currentThread], vNodes, vNodes, K_localArray[currentThread])
    end
    ### Multi-Thread Assembly ends here
    for thread ∈ 2:numOfThreads ### Adding all the thread vectors & matrices to 1st K_COO
        f[1] += f[thread]
        append!(K_COO[1].I, K_COO[thread].I)
        append!(K_COO[1].J, K_COO[thread].J)
        append!(K_COO[1].V, K_COO[thread].V)
        merge!(FeSpace, FeSpaceThreaded[thread]) ##Merge all threaded FeSpace together
    end
    return f[1], SparseArrays.SparseMatrixCSC(K_COO[1])
end


"""Responsible for assembling of FEM Scalars. The kind of things that you
may use this function for are, parameters important for the mathematical model.
Examples of such parameters are area of the element, it's center of garvity, moment of Interia etc.
parameters could be anything depending  on what you would like to send to the local matrix
assembler.
This function is useful if the local matrix makes changes to the parameters variable.

    v::Array{Float64,1} = RapidFEM.assembleScalar(parameters, volAttrib, FeSpace, mesh, RapidFEM.localScalar, problemDim, activeDimensions)
"""
function assembleScalar!(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh, localVectorFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1],
    varArgs...; reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::Vector where T

    numOfThreads::Int64 = Threads.nthreads()    #Total number of threads running
    noOfElements::Int64 = getNoOfElements(mesh, attribute)
    f::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, numOfThreads)
    f_localArray::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, numOfThreads)
    Threads.@threads for thread ∈ 1:numOfThreads
        f[thread] = zeros(noOfElements*problemDim)
        f_localArray[thread] = zeros(0)
    end
    FeSpaceThreaded::Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1} =
    Array{Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 1}(undef, numOfThreads)
    FeSpaceThreaded[1] = FeSpace #First thread can be the same as the original FeSpace
    Threads.@threads for thread ∈ 2:numOfThreads
        FeSpaceThreaded[thread] = deepcopy(FeSpace) #Others need to be deepcopies
    end
    RangeDict = createDimRange()
    dimRange::StepRange{Int64,Int64} = getRange(RangeDict, activeDimensions)
    Threads.@threads for elementNo ∈ 1:length(mesh.Elements[attribute])
        currentThread::Int64 = Threads.threadid()
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        if length(f_localArray[currentThread]) != problemDim
            f_localArray[currentThread] = zeros(problemDim)
        else
            fill!(f_localArray[currentThread],0.0)
        end
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, reduction = reduction, 
        elementFunction = elementFunction, quadrature = quadrature)

        localVectorFunc(f_localArray[currentThread], parameters, problemDim, element, elementNo, shapeFunction, coordArray, varArgs...)
        f[currentThread][problemDim*(elementNo-1)+1:problemDim*elementNo] += f_localArray[currentThread]
    end
    for thread ∈ 2:numOfThreads ### Adding all the thread vectors to 1st f
        f[1] += f[thread]
        merge!(FeSpace, FeSpaceThreaded[thread]) ##Merge all threaded FeSpace together
    end
    return f[1]
end

"""Responsible for assembling of FEM Scalars. The kind of things that you
may use this function for are, parameters important for the mathematical model.
Examples of such parameters are area of the element, it's center of garvity, moment of Interia etc.
parameters could be anything depending  on what you would like to send to the local matrix
assembler.

    v::Array{Float64,1} = RapidFEM.assembleScalar(parameters, volAttrib, FeSpace, mesh, RapidFEM.localScalar, problemDim, activeDimensions)
"""
function assembleScalar(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh, localVectorFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1], varArgs...; 
    reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::Vector where T

    return assembleScalar!(parameters, attribute, FeSpace,
    mesh, localVectorFunc, problemDim, activeDimensions, varArgs...;
    reduction = reduction, elementFunction = elementFunction, quadrature = quadrature)
end
