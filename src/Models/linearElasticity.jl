#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

function getTensorMapping()::Dict{Int64, Int64}
    mapDict::Dict{Int64, Int64} = Dict{Int64, Int64}()
    mapDict[11] = 1
    mapDict[22] = 2
    mapDict[33] = 3
    mapDict[12] = 4
    mapDict[23] = 5
    mapDict[13] = 6
    #mapDict[21] = 4
    #mapDict[32] = 5
    #mapDict[31] = 6
    return mapDict
end

function createVoigtElasticTensor(E::Float64, ν::Float64)::Array{Float64, 2}
    c::Float64 = E/((1+ν)*(1-2*ν))
    C::Array{Float64, 2} = zeros(6,6)
    C = [1-ν ν ν 0 0 0;
         ν 1-ν ν 0 0 0;
         ν ν 1-ν 0 0 0;
         0 0 0 (1-2*ν)/2 0 0;
         0 0 0 0 (1-2*ν)/2 0;
         0 0 0 0 0 (1-2*ν)/2]
    C = c*C
    return C
end

function getVoigtIndex(mapDict::Dict{Int64, Int64}, i::Int64, j::Int64)::Int64
    return mapDict[10*i+j]
end

function local_∇v_C_∇u!(K::Array{Float64,2},
    tensorMapN_ElasticTensor::Tuple{Dict{Int64, Int64}, Array{Float64, 2}},
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)
    mapDict::Dict{Int64, Int64} = tensorMapN_ElasticTensor[1]
    C::Array{Float64, 2} = tensorMapN_ElasticTensor[2]
    
    noOfIpPoints::Int64 = getNoOfElementIpPoints(shapeFunction)
    noOfNodes::Int64 = getNoOfElementNodes(shapeFunction)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
        ∂ϕ_∂x = get_∂ϕ_∂x(element, ∂x_∂ξ, shapeFunction, ipNo)
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
                                K[problemDim*(a-1)+i,problemDim*(b-1)+k] += c1*c2*∂ϕ_∂x[a,j]*C[ij,kl]*∂ϕ_∂x[b,l]*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+l] += c1*c2*∂ϕ_∂x[a,i]*C[ij,kl]*∂ϕ_∂x[b,k]*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+k] += c1*c2*∂ϕ_∂x[a,i]*C[ij,kl]*∂ϕ_∂x[b,l]*dΩ
                                K[problemDim*(a-1)+i,problemDim*(b-1)+l] += c1*c2*∂ϕ_∂x[a,j]*C[ij,kl]*∂ϕ_∂x[b,k]*dΩ
                            end
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function gaussianStress(tensorMapN_ElasticTensor::Tuple{Dict{Int64, Int64}, Array{Float64, 2}},
    solAtNodes::Array{Float64,1}, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)::Array{Array{Float64,1},1}
    mapDict::Dict{Int64, Int64} = tensorMapN_ElasticTensor[1]
    C::Array{Float64, 2} = tensorMapN_ElasticTensor[2]
    StressDim::Int64 = size(C,1)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    σ_g::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef, noOfIpPoints)
    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[1].∂ϕ_∂ξ)
    ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
    ∂ϕ_∂x::Array{Float64} = shapeFunction[1].∂ϕ_∂ξ*∂ξ_dx
    for ipNo::Int64 ∈ 1:noOfIpPoints
        σ_g[ipNo] = zeros(StressDim)
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx = ∂ξ_∂xFunc(∂x_∂ξ)
        ∂ϕ_∂x .= shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        for b::Int64 ∈ 1:noOfNodes
            for l::Int64 ∈ 1:problemDim
                for k::Int64 ∈ 1:l
                    kl::Int64 = RapidFEM.getVoigtIndex(mapDict, k, l)
                    c2::Float64 = (k==l) ? 0.5 : 1.0
                    for j::Int64 ∈ 1:problemDim
                        for i::Int64 ∈ 1:j
                            ij::Int64 = RapidFEM.getVoigtIndex(mapDict, i, j)
                            σ_g[ipNo][ij] += c2*C[ij,kl]*∂ϕ_∂x[b,l]*solAtNodes[problemDim*(b-1)+k]
                            σ_g[ipNo][ij] += c2*C[ij,kl]*∂ϕ_∂x[b,k]*solAtNodes[problemDim*(b-1)+l]
                        end
                    end
                end
            end
        end
    end
    return σ_g
end
