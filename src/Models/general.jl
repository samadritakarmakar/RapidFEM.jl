function local_∇v_∇u(parameterFunction::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})::Array{Float64,2}
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    K::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
        λ::Array{Float64, 1} = parameterFunction(x)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        for b ∈ 1:noOfNodes
            for a ∈ 1:noOfNodes
                for j ∈ 1:size(∂ϕ_∂x,2)
                    for i ∈ 1:problemDim
                        K[problemDim*(a-1)+i,problemDim*(b-1)+i] += λ[i]*∂ϕ_∂x[a,j]*∂ϕ_∂x[b,j]*dΩ
                        #println("K ", problemDim*(a-1)+i," ,", problemDim*(b-1)+i, " = ", λ[i]*∂ϕ_∂x[a,j]*∂ϕ_∂x[b,j]*dΩ)
                    end
                end
            end
        end
    end
    return K
end

function localSource(sourceFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    S::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        x::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        s::Array{Float64,1} = sourceFunc(x)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                S[problemDim*(a-1)+i] += ϕ[a]*s[i]*dΩ
            end
        end
    end
    return S
end

function localNeumann(neumannFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    #∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    Nm::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        #∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dS::Float64 = dSFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        x::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        nm::Array{Float64,1} = neumannFunc(x)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                Nm[problemDim*(a-1)+i] += ϕ[a]*nm[i]*dS
            end
        end
    end
    return Nm
end


function localScalar(scalarFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    S::Array{Float64,1} = zeros(problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        x::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        s::Array{Float64,1} = scalarFunc(x)
        for i ∈ 1:problemDim
            S[i] += s[i]*dΩ
        end
    end
    return S
end

function localScalarNeumann(scalarFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    #∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    S::Array{Float64,1} = zeros(problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        #∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dS::Float64 = dSFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        x::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        s::Array{Float64,1} = scalarFunc(x)
        for i ∈ 1:problemDim
            S[i] += s[i]*dS
        end
    end
    return S
end

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
    for ipNo ∈ 1:noOfIpPoint
        for node ∈ 1:numOfNodes
            f_n[fDim*(node-1)+1:fDim*node] += InvDist[ipNo][node]*f_g[ipNo]
        end
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

function getSolAtElement(sol::Array{Float64,1}, element::AbstractElement, problemDim::Int)::Array{Float64,1}
    vectorNodes::Array{Int64,1} = getVectorNodes(element, problemDim)
    solAtNodes::Array{Float64,1} = Array{Float64,1}(undef, length(vectorNodes))
    i::Int64 = 1
    for node ∈ vectorNodes
        solAtNodes[i] = sol[node]
        i +=1
    end
    return solAtNodes
end

function InvDistInterpolation(postProcessFunction::Function, sol::Array{Float64,1},
    parametersData::T,  FeSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}},
    mesh::Mesh,  attribute::Tuple{Int64, Int64}, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1], pow::Float64=12.0) where T
    f::Array{Float64,1} = []
    fDim::Int64 = 1;
    sumInvDistances::Vector = zeros(mesh.noOfNodes)
    RangeDict = RapidFEM.createDimRange()
    dimRange::StepRange{Int64,Int64} = RapidFEM.getRange(RangeDict, activeDimensions)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        solAtNodes::Array{Float64, 1} = getSolAtElement(sol, element, problemDim)
        coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
        coordArray::Array{Float64,2} = coordArrayTemp[dimRange,:]
        shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpace, element, mesh, lagrange)
        InvDists::Array{Array{Float64,1},1} = findInvDistances(coordArray, shapeFunction, pow)
        f_g::Array{Array{Float64,1},1} = postProcessFunction(parametersData,
        solAtNodes, problemDim, element, shapeFunction, coordArray)
        noOfElementNodes::Int64 = element.noOfElementNodes
        fDim = length(f_g[1])
        vectorFNodes::Array{Int64,1} = RapidFEM.getVectorNodes(element, fDim)
        if elementNo == 1
            fTemp::Array{Float64,1} = getWeighed_f(InvDists, f_g)
            f = zeros(fDim*mesh.noOfNodes)
            f[vectorFNodes] = fTemp
        else
            f[vectorFNodes] = getWeighed_f(InvDists, f_g)
        end
        nodes::Array{Int64,1} = RapidFEM.getNodes(element)
        sumInvDistances[nodes] = findSumInvDistances(InvDists)
    end
    finalizeInvInterpolation!(f,sumInvDistances)
    return f
end

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
