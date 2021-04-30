using RapidFEM, LinearAlgebra, Tensors

abstract type TrialTestFunction <:AbstractArray end

mutable struct TestFunction 
    mesh::RapidFEM.Mesh
    problemDim::Int64
    elementFunction::Function
    feSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}
    activeDimensions::Array{Int64, 1}
    currrentElementNo::Int64
    currentIpNo::Int64
    currentAttrib::Tuple{Int64, Int64}
    currentElement::RapidFEM.AbstractElement
    currentShapeFunction::Array{ShapeFunction, 1}
    function TestFunction(mesh::RapidFEM.Mesh, problemDim::Int64, elementFunction::Function, 
        feSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, activeDimensions::Array{Int64, 1} = [1,1,1])

        element = RapidFEM.PointElement(0, [0], [0], 0, 0)
        new(mesh, problemDim, elementFunction, feSpace, activeDimensions, 0, 0, (0,0), element, Array{ShapeFunction, 1}(undef, 0))
    end
end


mutable struct TrialFunction 
    testFunction::TestFunction
    lastSolution::Array{Float64}
    function TrialFunction(testFunction::TestFunction)
        new(testFunction, nonLinear, zeros(testFunction.problemDim*testFunction.mesh.noOfNodes), 0)
    end
end

function getNoOfCurrentElementNodes(v::TestFunction)
    return v.currentElement.noOfElementNodes
end

function getNoOfCurrentElementNodes(u::TrialFunction)
    return getNoOfCurrentElementNodes(u.testFunction)
end


function size(u_v::TrialTestFunction)
    return (getNoOfCurrentElementNodes(u_v)*u_v.problemDim,)
end

function getindex(v::TestFunction, nodeNo::Int)
    return  v.currentShapeFunction[v.currentIpNo].ϕ[nodeNo]*ones(v.problemDim)
end

function getindex(v::TestFunction, nodeNos...)
    returnVal = Array{Float64, 2}(v.problemDim, length(nodeNos))
    for nodeNo ∈ nodeNos
        returnVal[:, nodeNo] = v.currentShapeFunction[v.currentIpNo].ϕ[nodeNo]*ones(v.problemDim)
    end
    return  returnVal
end

function getindex(u::TrialFunction, nodeNo::Int)
    return  getindex(u.testFunction, nodeNo::Int)
end

function getindex(u::TrialFunction, nodeNo...)
    return  getindex(u.testFunction, nodeNo)
end

function setAttribute!(v::TestFunction, attribute::Tuple{Int64, Int64})
    v.currentAttrib = attribute
end

function setElement!(v::TestFunction, elementNo::Int64, reduction::Int64= 0)
    v.currentElementNo = elementNo
    v.currentElement = mesh.Elements[v.currentAttrib...][elementNo]
    v.currentShapeFunction = RapidFEM.feSpace!(v.feSpace, v.currentElement, v.mesh, reduction, v.elementFunction)
end

function setIpNo!(v::TestFunction, IpNo::Int64)
    v.currentIpNo = IpNo
end

abstract type TrialTestOperator end

mutable struct Grad <:TrialTestOperator
    u_v::TrialTestFunction
end






