using RapidFEM
function get_∂u_∂X!(∂u_∂X::Array{Float64, 2}, solAtNodes::Array{Float64, 1}, ∂ϕ_∂X::Array{Float64,2}, problemDim::Int64)
    fill!(∂u_∂X, 0.0)
    for a ∈ 1:size(∂ϕ_∂X, 1)
        for J ∈ 1:size(∂ϕ_∂X, 2)
            for i ∈ 1:problemDim
                ∂u_∂X[i,J] += ∂ϕ_∂X[a,J]*solAtNodes[problemDim*(a-1)+i]
            end
        end
    end
end

function local_δE_S_Vector!(f::Vector, hyperElasticData::T, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T


    model::HyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]
    lastSoln::Array{Float64, 1} = hyperElasticData[3]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)

    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    ∂u_∂X_array::Array{Float64, 2} = zeros(problemDim, problemDim)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X_array::Array{Float64} = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = LargeDefs.get1DTensor(∂ϕ_∂X_array')
        get_∂u_∂X!(∂u_∂X_array, solAtNodes, ∂ϕ_∂X_array, problemDim)
        ∂u_∂X = LargeDefs.get_∂u_∂X_Tensor(∂u_∂X_array)
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams)
        S_dot_Fᵀ = S⋅F'
        F_dot_S = F⋅S
        for a ∈ 1:noOfNodes
            ∂ϕ_∂X_a = ∂ϕ_∂X[a]
            for i ∈ 1:problemDim
                f[problemDim*(a-1)+i] += 0.5*(∂ϕ_∂X_a ⊗ F[i,:] + F[i,:] ⊗ ∂ϕ_∂X_a) ⊡ S*dΩ
            end
        end
    end
    return nothing
end

function local_δE_Cᵀ_ΔE!(𝕂::Array{Float64,2}, hyperElasticData::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T


    model::HyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]
    lastSoln::Array{Float64, 1} = hyperElasticData[3]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    ∂u_∂X_array::Array{Float64, 2} = zeros(problemDim, problemDim)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = RapidFEM.get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        #ϕ = shapeFunction[ipNo].ϕ
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X_array = RapidFEM.get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = LargeDefs.get1DTensor(∂ϕ_∂X_array')
        get_∂u_∂X!(∂u_∂X_array, solAtNodes, ∂ϕ_∂X_array, problemDim)
        ∂u_∂X = LargeDefs.get_∂u_∂X_Tensor(∂u_∂X_array)
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S = model.secondPiolaStress(E, modelParams)
        ℂ = model.materialTangentTensor(E, modelParams)
        for b::Int64 ∈ 1:noOfNodes
            ∂ϕ_∂X_b = ∂ϕ_∂X[b]
            for a::Int64 ∈ 1:noOfNodes
                ∂ϕ_∂X_a = ∂ϕ_∂X[a]
                for j ∈ 1:problemDim
                    𝕂[problemDim*(a-1)+j, problemDim*(b-1)+j] += (∂ϕ_∂X_a⋅S⋅∂ϕ_∂X_b)*dΩ
                    F_j = F[j,:]
                    for i ∈ 1:problemDim
                        #𝕂[problemDim*(a-1)+i, problemDim*(b-1)+j] +=
                        #0.25*(∂ϕ_∂X[a] ⊗ F[i,:] + F[i,:] ⊗ ∂ϕ_∂X[a]) ⊡ ℂ ⊡
                        #(∂ϕ_∂X[b] ⊗ F_j + F_j ⊗ ∂ϕ_∂X[b])*dΩ
                        F_i = F[i,:]
                        𝕂[problemDim*(a-1)+i, problemDim*(b-1)+j] +=
                        0.25*(∂ϕ_∂X_a ⊗ F_i + F_i ⊗ ∂ϕ_∂X_a) ⊡ ℂ ⊡
                        (∂ϕ_∂X_b ⊗ F_j + F_j ⊗ ∂ϕ_∂X_b)*dΩ
                    end
                end
            end
        end
    end
    #SmallStrainPlastic.updateStateDict!(stateDictBuffer, stateDictBufferCopy)
    return nothing
end


function localReferenceSource!(S::Vector, hyperElasticData::T, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    sourceFunc = hyperElasticData[1]
    lastSoln::Array{Float64, 1} = hyperElasticData[2]

    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    #S::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ::Float64 = RapidFEM.get_dΩ(element, ∂X_∂ξ, shapeFunction, ipNo)
        ϕ::Array{Float64,1} = RapidFEM.get_ϕ(shapeFunction, ipNo)
        X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        s::Array{Float64,1} = sourceFunc(X; kwargs4function...)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                S[problemDim*(a-1)+i] += ϕ[a]*s[i]*dΩ
            end
        end
    end
    return nothing
end


function localReferenceNeumann!(Nm::Vector, hyperElasticData::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    neumannFunc = hyperElasticData[1]
    lastSoln::Array{Float64, 1} = hyperElasticData[2]

    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    #Nm::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} =  get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dS::Float64 = RapidFEM.get_dS(element, ∂X_∂ξ,  shapeFunction, ipNo)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        nm::Array{Float64,1} = neumannFunc(X; kwargs4function...)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                Nm[problemDim*(a-1)+i] += ϕ[a]*nm[i]*dS
            end
        end
    end
    return nothing
end

function gaussianSecondPiolaStress(hyperElasticData::T, solAtNodes::Array{Float64,1}, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    model::HyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]

    
    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    ∂u_∂X_array::Array{Float64, 2} = zeros(problemDim, problemDim)
    S_g::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, noOfIpPoints)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} =  get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X_array::Array{Float64} = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂X = LargeDefs.get1DTensor(∂ϕ_∂X_array')
        get_∂u_∂X!(∂u_∂X_array, solAtNodes, ∂ϕ_∂X_array, problemDim)
        ∂u_∂X = LargeDefs.get_∂u_∂X_Tensor(∂u_∂X_array)
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S_tensor = model.secondPiolaStress(E, modelParams)
        S_g[ipNo] = zeros(problemDim^2)
        k = 1
        for i ∈ 1:problemDim
            for j ∈ 1:problemDim
                S_g[ipNo][k] = S_tensor[i,j]
                k +=1
            end
        end
    end
    return S_g
end

function gaussianGreenStrain(hyperElasticData::T, solAtNodes::Array{Float64,1}, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    model::HyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]

    
    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    ∂u_∂X_array::Array{Float64, 2} = zeros(problemDim, problemDim)
    E_g::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, noOfIpPoints)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} =  get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        #ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X_array::Array{Float64} = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        get_∂u_∂X!(∂u_∂X_array, solAtNodes, ∂ϕ_∂X_array, problemDim)
        ∂u_∂X = LargeDefs.get_∂u_∂X_Tensor(∂u_∂X_array)
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E_tensor = LargeDefs.getGreenLagrangeStrain(F)
        E_g[ipNo] = zeros(problemDim^2)
        k = 1
        for i ∈ 1:problemDim
            for j ∈ 1:problemDim
                E_g[ipNo][k] = E_tensor[i,j]
                k +=1
            end
        end
    end
    return E_g
end

function gaussianDeformationGrad(hyperElasticData::T, solAtNodes::Array{Float64,1}, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    model::HyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]

    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    ∂u_∂X_array::Array{Float64, 2} = zeros(problemDim, problemDim)
    F_g::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, noOfIpPoints)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} =  get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X_array::Array{Float64} = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        get_∂u_∂X!(∂u_∂X_array, solAtNodes, ∂ϕ_∂X_array, problemDim)
        ∂u_∂X = LargeDefs.get_∂u_∂X_Tensor(∂u_∂X_array)
        F_tensor = LargeDefs.getDeformationGradient(∂u_∂X)
        F_g[ipNo] = zeros(problemDim^2)
        k = 1
        for i ∈ 1:problemDim
            for j ∈ 1:problemDim
                F_g[ipNo][k] = F_tensor[i,j]
                k +=1
            end
        end
    end
    return F_g
end

function gaussianCauchyStress(hyperElasticData::T, solAtNodes::Array{Float64,1}, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    model::HyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]

    
    noOfIpPoints::Int64 = RapidFEM.getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = RapidFEM.getNoOfElementNodes(shapeFunction)
    ∂u_∂X_array::Array{Float64, 2} = zeros(problemDim, problemDim)
    σ_g::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, noOfIpPoints)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} =  get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X_array::Array{Float64} = get_∂ϕ_∂x(element, ∂X_∂ξ, shapeFunction, ipNo)
        get_∂u_∂X!(∂u_∂X_array, solAtNodes, ∂ϕ_∂X_array, problemDim)
        ∂u_∂X = LargeDefs.get_∂u_∂X_Tensor(∂u_∂X_array)
        F = LargeDefs.getDeformationGradient(∂u_∂X)
        E = LargeDefs.getGreenLagrangeStrain(F)
        S_tensor = model.secondPiolaStress(E, modelParams)
        σ_tensor = 1.0/det(F)*F⋅S_tensor⋅F'
        σ_g[ipNo] = zeros(problemDim^2)
        k = 1
        for i ∈ 1:problemDim
            for j ∈ 1:problemDim
                σ_g[ipNo][k] = σ_tensor[i,j]
                k +=1
            end
        end
    end
    return σ_g
end
