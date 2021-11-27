#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

function findInvDistances(coordArray::Array{Float64,2}, shapeFunction::Array{ShapeFunction,1}, pow::Float64)::Array{Array{Float64,1},1}
    noOfIpPoints::Int64 = length(shapeFunction)
    InvDists::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, noOfIpPoints)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
        InvDists[ipNo] = Array{Float64,1}(undef, size(coordArray,2))
        for node ∈ 1:size(coordArray,2)
            dist::Float64 = norm(coordArray[:,node] - x)
            InvDists[ipNo][node] = 1 ./ dist^pow
        end
    end
    return InvDists
end



function getWeighed_f(InvDist::Array{Array{Float64,1},1}, f_g::Array{Array{Float64,1},1})::Array{Float64,1}
    numOfNodes::Int64 = length(InvDist[1])
    fDim::Int64 = length(f_g[1])
    noOfIpPoint::Int64 = length(f_g)
    #println(numOfNodes, " ", fDim, " ", noOfIpPoint)
    f_n::Array{Float64,1} = zeros(numOfNodes*fDim)
    ipNo = 1
    InfFlag = false
    while (ipNo <= noOfIpPoint && InfFlag == false)
        for node ∈ 1:numOfNodes
            if InvDist[ipNo][node] == Inf
                f_n[fDim*(node-1)+1:fDim*node]  = f_g[ipNo]
                InfFlag = true
                break
            else
                f_n[fDim*(node-1)+1:fDim*node] += InvDist[ipNo][node]*f_g[ipNo]
            end
        end
        ipNo += 1
    end
    return f_n
end

function findSumInvDistances(InvDist::Array{Array{Float64,1},1})::Array{Float64,1}
    numOfNodes::Int64 = length(InvDist[1])
    noOfIpPoint::Int64 = length(InvDist)
    sumInvDist::Array{Float64,1} = zeros(numOfNodes)
    for ipNo ∈ 1:noOfIpPoint
        for node ∈ 1:numOfNodes
            sumInvDist[node] += InvDist[ipNo][node]
        end
    end
    return sumInvDist
end

function getFDim(f::Vector,element::AbstractElement)::Int64
    numOfNodes::Int64 = element.noOfElementNodes
    #println(length(f))
    fDim::Int64 = length(f)/numOfNodes
    return fDim
end

function finalizeInvInterpolation!(f::Vector,sumInvDistances::Vector)
    numOfNodes::Int64 = length(sumInvDistances)
    fDim::Int64 = length(f)/numOfNodes
    for node ∈ 1:numOfNodes
        f[fDim*(node-1)+1:fDim*node] = f[fDim*(node-1)+1:fDim*node]/sumInvDistances[node]
    end
end

"""Extracts the solution available at a particular Element for a certain problem dimension.

    solAtNodes::Array{Float64,1} = getSolAtElement(sol, element, problemDim)
"""
function getSolAtElement(sol::Array{Float64,1}, element::AbstractElement, problemDim::Int)::Array{Float64,1}
    vectorNodes::Array{Int64,1} = getVectorNodes(element, problemDim)
    #sort!(vectorNodes)
    solAtNodes::Array{Float64,1} = Array{Float64,1}(undef, length(vectorNodes))
    i::Int64 = 1
    for node ∈ vectorNodes
        solAtNodes[i] = sol[node]
        i +=1
    end
    return solAtNodes
end


"""This function is used to for post processing, meaning for obtaining additional data from the solution
obtained. One such use is to find the stress from the generated displacement data from the solution.

See src/Examples/LinearElastic2Material.jl example.

        σTemp::Array{Float64,1} = RapidFEM.InvDistInterpolation([RapidFEM.gaussianStress, RapidFEM.gaussianStress], x, [(tensorMap, C1), (tensorMap, C1)],  FeSpace, mesh,  [volAttrib1, volAttrib2], problemDim, activeDimensions)
"""
function InvDistInterpolation(postProcessFunctionArray::Array{func, 1}, sol::Array{Float64,1},
    parametersDataArray::Array{T, 1},  FeSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}},
    mesh::Mesh,  attributeArray::Array{Tuple{Int64, Int64},1}, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1], varArgs...; pow::Float64=14.0, reduction::Int64 = 0, elementFunction = lagrange, quadrature = gauss) where {func, T}
    @assert length(postProcessFunctionArray)==length(attributeArray) "Length of postProcessFunctionArray should be equal to attributeArray."
    @assert length(parametersDataArray)==length(attributeArray) "Length of parametersDataArray should be equal to attributeArray."
    f::Array{Float64,1} = []
    fDim::Int64 = 1;
    sumInvDistances::Vector = zeros(mesh.noOfNodes)
    RangeDict = RapidFEM.createDimRange()
    dimRange::StepRange{Int64,Int64} = RapidFEM.getRange(RangeDict, activeDimensions)
    for attributeNo ∈ 1:length(attributeArray)
        attribute::Tuple{Int64, Int64} = attributeArray[attributeNo]
        postProcessFunction::func = postProcessFunctionArray[attributeNo]
        parametersData::T = parametersDataArray[attributeNo]
        for elementNo ∈ 1:length(mesh.Elements[attribute])
            element::AbstractElement = mesh.Elements[attribute][elementNo]
            solAtNodes::Array{Float64, 1} = getSolAtElement(sol, element, problemDim)
            coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
            coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
            shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, reduction = reduction, 
            elementFunction = elementFunction, quadrature = quadrature)
            InvDists::Array{Array{Float64,1},1} = findInvDistances(coordArray, shapeFunction, pow)
            f_g::Array{Array{Float64,1},1} = postProcessFunction(parametersData,
            solAtNodes, problemDim, element, elementNo, shapeFunction, coordArray; varArgs...)
            noOfElementNodes::Int64 = element.noOfElementNodes
            fDim = length(f_g[1])
            vectorFNodes::Array{Int64,1} = RapidFEM.getVectorNodes(element, fDim)
            if (elementNo == 1 && attributeNo ==1)
                fTemp::Array{Float64,1} = getWeighed_f(InvDists, f_g)
                f = zeros(fDim*mesh.noOfNodes)
                f[vectorFNodes] += fTemp
            else
                f[vectorFNodes] += getWeighed_f(InvDists, f_g)
            end
            nodes::Array{Int64,1} = RapidFEM.getNodes(element)
            sumInvDistances[nodes] += findSumInvDistances(InvDists)
        end
    end
    for node ∈ 1:mesh.noOfNodes #Check no Inv Distance is a zero, if so replace with 1.0
        sumInvDistances[node] = sumInvDistances[node]==Inf ? 1.0 : sumInvDistances[node]
    end
    finalizeInvInterpolation!(f,sumInvDistances)
    return f
end

"""In continnumm mechanics many terms are symmetric in nature. Such as Cauchy Stresses
These terms are calculated in their voigt notations. This function transforms those
terms back to tensor notations for representation in VTK files. As of now only 3x3 2nd order
tensors are supported.

See src/Examples/LinearElastic2Material.jl example.

    σ::Array{Float64,1} = RapidFEM.voigtToTensor(σTemp, mesh)
"""
function voigtToTensor(sol::Array{Float64,1}, mesh::Mesh)::Array{Float64,1}
    noOfNodes::Int64 = mesh.noOfNodes
    I::Array{Int64,1} = collect(1:9)
    J::Array{Int64,1} = [1, 4, 6, 4, 2, 5, 6, 5, 3]
    V::Array{Float64,1} = ones(9)
    P::SparseMatrixCSC = sparse(I,J,V)
    newSol::Array{Float64,1} = zeros(noOfNodes*9)
    for node ∈ 1:noOfNodes
        newSol[9*(node-1)+1:9*node] = P*sol[6*(node-1)+1:6*node]
    end
    return newSol
end
