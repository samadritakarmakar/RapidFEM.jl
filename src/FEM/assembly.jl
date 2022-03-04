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

function getProblemType(meshTuple::Tuple{Mesh, Mesh}, problemDimTuple::Tuple{Int64, Int64})
    if meshTuple[1] == meshTuple[2] && problemDimTuple[1] == problemDimTuple[2]
        return :General
    elseif meshTuple[1] == meshTuple[2]
        return :Coupled
    else
        return :Mixed
    end  
end

function getElements(problemType::Symbol, attribute::Tuple{Int64, Int64}, meshTuple::Tuple{Mesh, Mesh}, elementNo::Int64)
    
    element_2::AbstractElement =  meshTuple[2].Elements[attribute][elementNo]
    if problemType == :Mixed
        element_1::AbstractElement = meshTuple[1].Elements[attribute][elementNo]
        return element_1, element_2
    end
    return element_2, element_2
end

function getReductionFix(problemType::Symbol, element_1::AbstractElement, element_2::AbstractElement)

    if problemType != :Mixed
        return 0, 0
    end
    #By default element_2 is the set order. If element_1 has higher order, then that is used.
    element = element_2
    order = element_2.order
    if element_1.order > element_2.order
        element = element_1
        order = element_1.order          
    end
    #Fix order based on: order_max = order_i - reduction_fix_i; reduction_fix_i is negative since a negative reduction means an increase in order.
    reduction_fix_1 = element_1.order - order
    reduction_fix_2 = element_2.order - order
    return reduction_fix_1, reduction_fix_2
end

"""Responsible for assembling of FEM Matrices. paramerterFunction could
be anything depending on what you would like to send to the local matrix
assembler. This function is to be used if the Matrix shall be Mixed and Coupled.
This function is useful if the local matrix makes changes to the parameters variable.

    K::SparseMatrixCSC = assembleMatrix!(parameters, attribute, FeSpace, (mesh1, mesh2), localMatrixFunc, (problemDim1, problemDim2),
    activeDimensions, varArgs...; reduction = reduction, elementFunction=lagrange, quadrature= gauss)
"""
function assembleMatrix!(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    meshTuple::Tuple{Mesh, Mesh}, localMatrixFunc::Function, problemDimTuple::Tuple{Int64, Int64},
    activeDimensions::Array{Int64,1}=[1, 1, 1], varArgs...; 
    reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::SparseMatrixCSC where T   

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

    problemType = getProblemType(meshTuple, problemDimTuple)

    #######This makes sure the full size of K_COO is same as the size of the filled K matrix
    vNodes2_1::Array{Int64,1} = [meshTuple[1].noOfNodes*problemDimTuple[1]]
    vNodes2_2::Array{Int64,1} = [meshTuple[2].noOfNodes*problemDimTuple[2]]
    K_local2::Array{Float64,2} = zeros(1,1)
    ####
    Threads.@threads for thread ∈ 1:numOfThreads #Setting the matrices size for each thread
        FEMSparse.assemble_local_matrix!(K_COO[thread], vNodes2_1, vNodes2_2, K_local2)
    end
    ### Multi-Thread Assembly starts here
    Threads.@threads for elementNo ∈ 1:length(meshTuple[2].Elements[attribute])
        currentThread::Int64 = Threads.threadid()
        #Element Initializations
        element_1, element_2 = getElements(problemType, attribute, meshTuple, elementNo)
        #if the type of element changes then reallocate the local matrix else replace with zeros
        if size(K_localArray[currentThread]) != (problemDimTuple[1]*element_1.noOfElementNodes, problemDimTuple[2]*element_2.noOfElementNodes)
            K_localArray[currentThread] = zeros(problemDimTuple[1]*element_1.noOfElementNodes, problemDimTuple[2]*element_2.noOfElementNodes)
        else
            fill!(K_localArray[currentThread],0.0)
        end
        reduction_fix_1, reduction_fix_2 = getReductionFix(problemType, element_1, element_2)
        shapeFunction_2::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element_2, meshTuple[2], reduction = (reduction + reduction_fix_2), 
        elementFunction = elementFunction, quadrature = quadrature)
        #Arguments for localMatrixFunc
        coordArrayTemp2::Array{Float64,2} = getCoordArray(meshTuple[2], element_2)
        coordArray2::Array{Float64,2} = coordArrayTemp2[dimRange,:]
        if problemType == :General
            localMatrixFunc(K_localArray[currentThread], parameters, problemDimTuple[2], element_2, elementNo, shapeFunction_2, coordArray2, varArgs...)
        elseif problemType == :Coupled
            localMatrixFunc(K_localArray[currentThread], parameters, problemDimTuple[1], problemDimTuple[2], element_2, elementNo, shapeFunction_2, coordArray2, varArgs...)
        elseif problemType == :Mixed
            shapeFunction_1::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element_1, meshTuple[1], reduction = (reduction + reduction_fix_1), 
                elementFunction = elementFunction, quadrature = quadrature)
            coordArrayTemp1::Array{Float64,2} = getCoordArray(meshTuple[1], element_1)
            coordArray1::Array{Float64,2} = coordArrayTemp1[dimRange,:]
            localMatrixFunc(K_localArray[currentThread], parameters, problemDimTuple[1], problemDimTuple[2], element_1, element_2, elementNo, shapeFunction_1, shapeFunction_2, coordArray1, coordArray2, varArgs...)
        else
            error("Unknown Problem Type!")
        end
        vNodes_1::Array{Int64} = getVectorNodes(element_1, problemDimTuple[1])
        vNodes_2::Array{Int64} = getVectorNodes(element_2, problemDimTuple[2])
        FEMSparse.assemble_local_matrix!(K_COO[currentThread], vNodes_1, vNodes_2, K_localArray[currentThread])
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
function assembleMatrix!(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh, localMatrixFunc::Function, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1], varArgs...; 
    reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::SparseMatrixCSC where T

    return assembleMatrix!(parameters, attribute, FeSpace, (mesh, mesh), localMatrixFunc, (problemDim, problemDim),
    activeDimensions, varArgs...; reduction = reduction, elementFunction=elementFunction, quadrature= quadrature)
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

    return assembleMatrix!(parameters, attribute, FeSpace, mesh,  localMatrixFunc, problemDim,
    activeDimensions, varArgs...; reduction = reduction, elementFunction = elementFunction, quadrature = quadrature)
end

"""Responsible for assembling of FEM Matrices. paramerterFunction could
be anything depending on what you would like to send to the local matrix
assembler. This function is to be used if the Matrix shall be Coupled.
This function is useful if the local matrix makes changes to the parameters variable.

    K::SparseMatrixCSC = assembleMatrix!(parameters, volAttrib, FeSpace, mesh, RapidFEM.local_∇v_λ_∇u!, problemDim, activeDimensions)
"""
function assembleMatrix!(parameters::T, attribute::Tuple{Int64, Int64},
    FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh, localMatrixFunc::Function, problemDimTuple::Tuple{Int64, Int64},
    activeDimensions::Array{Int64,1}=[1, 1, 1], varArgs...; 
    reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::SparseMatrixCSC where T

    return assembleMatrix!(parameters, attribute, FeSpace, (mesh, mesh), localMatrixFunc, problemDimTuple,
    activeDimensions, varArgs...; reduction = reduction, elementFunction=elementFunction, quadrature= quadrature)
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
