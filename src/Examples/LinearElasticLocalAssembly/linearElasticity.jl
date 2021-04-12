#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

function createElasticTensor(E::Float64, ν::Float64)
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    C = λ*one(SymmetricTensor{2,3, Float64})⊗ one(SymmetricTensor{2,3, Float64})
    C +=2*μ*one(SymmetricTensor{4,3, Float64})
    return C
end


function local_∇v_C_∇u_Tensor!(K::Array{Float64,2},
    tensorMapN_ElasticTensor::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T
    C = tensorMapN_ElasticTensor[1]
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[1].∂ϕ_∂ξ)
    ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
    x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[1].ϕ)
    dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[1].ipData)
    ∂ϕ_∂x::Array{Float64} = shapeFunction[1].∂ϕ_∂ξ*∂ξ_dx
    #K::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx = ∂ξ_∂xFunc(∂x_∂ξ)
        x = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
        dΩ = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ∂ϕ_∂x .= shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        for b::Int64 ∈ 1:noOfNodes
            for a::Int64 ∈ 1:noOfNodes
                for l::Int64 ∈ 1:problemDim
                    for k::Int64 ∈ 1:l
                        for j::Int64 ∈ 1:problemDim
                            for i::Int64 ∈ 1:j
                                K[problemDim*(a-1)+i,problemDim*(b-1)+k] += 0.25*∂ϕ_∂x[a,j]*C[i,j,k,l]*∂ϕ_∂x[b,l]*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+l] += 0.25*∂ϕ_∂x[a,i]*C[i,j,k,l]*∂ϕ_∂x[b,k]*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+k] += 0.25*∂ϕ_∂x[a,i]*C[i,j,k,l]*∂ϕ_∂x[b,l]*dΩ
                                K[problemDim*(a-1)+i,problemDim*(b-1)+l] += 0.25*∂ϕ_∂x[a,j]*C[i,j,k,l]*∂ϕ_∂x[b,k]*dΩ
                            end
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function gaussianStress(tensorMapN_ElasticTensor::T,
    solAtNodes::Array{Float64,1}, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T
    C = tensorMapN_ElasticTensor[1]
    StressDim::Int64 = size(C,1)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    σ_g::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, noOfIpPoints)
    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[1].∂ϕ_∂ξ)
    ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
    ∂ϕ_∂x::Array{Float64} = shapeFunction[1].∂ϕ_∂ξ*∂ξ_dx
    for ipNo::Int64 ∈ 1:noOfIpPoints
        σ_g[ipNo] = zeros(StressDim^2)
        for b::Int64 ∈ 1:noOfNodes
            for l::Int64 ∈ 1:problemDim
                for k::Int64 ∈ 1:problemDim
                    ij = 1
                    for j::Int64 ∈ 1:problemDim
                        for i::Int64 ∈ 1:problemDim
                            σ_g[ipNo][ij] += 0.5*C[i,j,k,l]*∂ϕ_∂x[b,l]*solAtNodes[problemDim*(b-1)+k]
                            σ_g[ipNo][ij] += 0.5*C[i,j,k,l]*∂ϕ_∂x[b,k]*solAtNodes[problemDim*(b-1)+l]
                            ij +=1
                        end
                    end
                end
            end
        end
    end
    return σ_g
end
