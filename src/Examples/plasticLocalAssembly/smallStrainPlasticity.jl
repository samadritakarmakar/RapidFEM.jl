#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

function findStrain!(mapDict::Dict{Int64, Int64}, ϵ::Array{Float64, 1}, ∂ϕ_∂x::Array{Float64,2},
    solAtNodes::Array{Float64, 1}, problemDim::Int64)

    fill!(ϵ, 0.0)
    for a ∈ 1:size(∂ϕ_∂x, 1)
        for j ∈ 1:problemDim
            for i ∈ 1:j
                ij::Int64 = getVoigtIndex(mapDict, i, j)
                c1::Float64 = (i==j) ? 0.5 : 1.0
                ϵ[ij] += c1*(∂ϕ_∂x[a,j]*solAtNodes[problemDim*(a-1)+i]+∂ϕ_∂x[a,i]*solAtNodes[problemDim*(a-1)+j])
            end
        end
    end
end

function local_∇v_Cᵀ_∇u!(K::Array{Float64,2}, tensorMap_N_PlasticData::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T

    mapDict::Dict{Int64, Int64} = tensorMap_N_PlasticData[1]
    C::Array{Float64,2} = tensorMap_N_PlasticData[2]
    model::PlasticModel = tensorMap_N_PlasticData[3]
    modelParams::ModelParams  = tensorMap_N_PlasticData[4]
    stateDict  = tensorMap_N_PlasticData[5]
    stateDictBuffer  = tensorMap_N_PlasticData[6]
    stateDictBufferCopy = deepcopy(stateDictBuffer)
    lastSoln::Array{Float64, 1} = tensorMap_N_PlasticData[7]
    plasticVars::PlasticVars = SmallStrainPlastic.initPlasticVars(model)
    plasticVars.C = C

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    ϵ::Array{Float64, 1} = zeros(model.ϵSize)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        findStrain!(mapDict, ϵ, ∂ϕ_∂x,  solAtNodes, problemDim)
        plasticVars.ϵ = deepcopy(ϵ)
        getState!(plasticVars.ϵᵖ, plasticVars.α, stateDictBufferCopy, elementNo, ipNo)
        plasticVars.Cᵀ = SmallStrainPlastic.findNumerical_Cᵀ(plasticVars, model,
        modelParams, stateDict,  elementNo, ipNo)
        #SmallStrainPlastic.checkPlasticState!(plasticVars, model,
        #modelParams, stateDict, stateDictBufferCopy,  elementNo, ipNo)
        #println("plasticVars.ϵᵖ at ∇v_Cᵀ_∇u= ", plasticVars.ϵᵖ)
        #println("plasticVars.Cᵀ = ", plasticVars.Cᵀ)
        for b::Int64 ∈ 1:noOfNodes
            for a::Int64 ∈ 1:noOfNodes
                for l::Int64 ∈ 1:problemDim
                    for k::Int64 ∈ 1:l
                        kl::Int64 = getVoigtIndex(mapDict, k, l)
                        c2::Float64 = (k==l) ? 0.5 : 1.0
                        for j::Int64 ∈ 1:problemDim
                            for i::Int64 ∈ 1:j
                                ij::Int64 = getVoigtIndex(mapDict, i, j)
                                c1::Float64 = (i==j) ? 0.5 : 1.0
                                K[problemDim*(a-1)+i,problemDim*(b-1)+k] += c1*c2*∂ϕ_∂x[a,j]*plasticVars.Cᵀ[ij,kl]*∂ϕ_∂x[b,l]*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+l] += c1*c2*∂ϕ_∂x[a,i]*plasticVars.Cᵀ[ij,kl]*∂ϕ_∂x[b,k]*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+k] += c1*c2*∂ϕ_∂x[a,i]*plasticVars.Cᵀ[ij,kl]*∂ϕ_∂x[b,l]*dΩ
                                K[problemDim*(a-1)+i,problemDim*(b-1)+l] += c1*c2*∂ϕ_∂x[a,j]*plasticVars.Cᵀ[ij,kl]*∂ϕ_∂x[b,k]*dΩ
                            end
                        end
                    end
                end
            end
        end
    end
    SmallStrainPlastic.updateStateDict!(stateDictBuffer, stateDictBufferCopy)
    return nothing
end

function local_∇v_σ_Vector!(f::Vector, tensorMap_N_PlasticData::T, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T

    mapDict::Dict{Int64, Int64} = tensorMap_N_PlasticData[1]
    C::Array{Float64,2} = tensorMap_N_PlasticData[2]
    model::PlasticModel = tensorMap_N_PlasticData[3]
    modelParams::ModelParams  = tensorMap_N_PlasticData[4]
    stateDict  = tensorMap_N_PlasticData[5]
    stateDictBuffer  = tensorMap_N_PlasticData[6]
    stateDictBufferCopy = deepcopy(stateDictBuffer)
    lastSoln::Array{Float64, 1} = tensorMap_N_PlasticData[7]
    plasticVars::PlasticVars = SmallStrainPlastic.initPlasticVars(model)
    plasticVars.C = C

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    ϵ::Array{Float64, 1} = zeros(model.ϵSize)
    for ipNo ∈ 1:noOfIpPoints
        #ϵ::Array{Float64, 1} = zeros(model.ϵSize)
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        x::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        findStrain!(mapDict, ϵ, ∂ϕ_∂x,  solAtNodes, problemDim)
        plasticVars.ϵ = deepcopy(ϵ)
        getState!(plasticVars.ϵᵖ, plasticVars.α, stateDictBufferCopy, elementNo, ipNo)
        SmallStrainPlastic.checkPlasticState!(plasticVars, model,
        modelParams, stateDict, stateDictBufferCopy,  elementNo, ipNo)
        #println("plasticVars.ϵᵖ ∇v_σ= \n", plasticVars.ϵᵖ)
        for a ∈ 1:noOfNodes
            for j::Int64 ∈ 1:problemDim
                for i::Int64 ∈ 1:j
                    ij::Int64 = getVoigtIndex(mapDict, i, j)
                    c1::Float64 = (i==j) ? 0.5 : 1.0
                    f[problemDim*(a-1)+i] += c1*∂ϕ_∂x[a,j]*plasticVars.σ_voigt[ij]*dΩ
                    f[problemDim*(a-1)+j] += c1*∂ϕ_∂x[a,i]*plasticVars.σ_voigt[ij]*dΩ
                end
            end
        end
    end
    SmallStrainPlastic.updateStateDict!(stateDictBuffer, stateDictBufferCopy)
    return nothing
end
#=
j2Model = SmallStrainPlastic.j2Model

function initParams_j2(σ_y::Float64, params_H::Float64)
    return SmallStrainPlastic.initParams_j2(σ_y, params_H)
end

function updateStateDict4rmBuffer()
    SmallStrainPlastic.updateStateDict4rmBuffer()
    return nothing
end
=#
