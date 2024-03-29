#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

function local_∇v_λ_∇u!(K::Array{Float64,2}, parameters::Function,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)

    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #K[:,:] = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[1].∂ϕ_∂ξ)
    ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
    x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[1].ϕ)
    λ::Array{Float64, 1} = parameters(x; kwargs4function...)
    dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[1].ipData)
    ∂ϕ_∂x::Array{Float64} = shapeFunction[1].∂ϕ_∂ξ*∂ξ_dx
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx = ∂ξ_∂xFunc(∂x_∂ξ)
        x = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
        λ = parameters(x; kwargs4function...)
        dΩ = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ∂ϕ_∂x .= shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
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
    return nothing
end

function local_v_ρ_u!(M::Array{Float64,2}, parameters::Function,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)

    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[1].∂ϕ_∂ξ)
    ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
    ϕ::Array{Float64} = shapeFunction[1].ϕ
    x::Array{Float64, 1} = getInterpolated_x(coordArray, ϕ)
    ρ::Array{Float64, 1} = parameters(x; kwargs4function...)
    dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[1].ipData)
    #M::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx = ∂ξ_∂xFunc(∂x_∂ξ)
        ϕ = shapeFunction[ipNo].ϕ
        x = getInterpolated_x(coordArray, ϕ)
        ρ = parameters(x; kwargs4function...)
        dΩ = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        for b ∈ 1:noOfNodes
            for a ∈ 1:noOfNodes
                for i ∈ 1:problemDim
                    M[problemDim*(a-1)+i,problemDim*(b-1)+i] += ρ[i]*ϕ[a]*ϕ[b]*dΩ
                    #println("K ", problemDim*(a-1)+i," ,", problemDim*(b-1)+i, " = ", λ[i]*∂ϕ_∂x[a,j]*∂ϕ_∂x[b,j]*dΩ)
                end
            end
        end
    end
    return nothing
end

function localBoundary_v_ρ_u!(M::Array{Float64,2}, parameters::Function,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)

    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
    ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
    ϕ::Array{Float64} = shapeFunction[1].ϕ
    x::Array{Float64, 1} = getInterpolated_x(coordArray, ϕ)
    ρ::Array{Float64, 1} = parameters(x; kwargs4function...)
    dS::Float64 = dSFunc(∂x_∂ξ, shapeFunction[1].ipData)
    #M::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx = ∂ξ_∂xFunc(∂x_∂ξ)
        ϕ = shapeFunction[ipNo].ϕ
        x = getInterpolated_x(coordArray, ϕ)
        ρ = parameters(x; kwargs4function...)
        dS = dSFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        for b ∈ 1:noOfNodes
            for a ∈ 1:noOfNodes
                for i ∈ 1:problemDim
                    M[problemDim*(a-1)+i,problemDim*(b-1)+i] += ρ[i]*ϕ[a]*ϕ[b]*dS
                    #println("K ", problemDim*(a-1)+i," ,", problemDim*(b-1)+i, " = ", λ[i]*∂ϕ_∂x[a,j]*∂ϕ_∂x[b,j]*dΩ)
                end
            end
        end
    end
    return nothing
end


function localSource!(S::Vector, sourceFunc::Function, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)

    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        #∂ξ_dx = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
        ϕ = get_ϕ(shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, ϕ)
        s = sourceFunc(x; kwargs4function...)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                S[problemDim*(a-1)+i] += ϕ[a]*s[i]*dΩ
            end
        end
    end
    return nothing
end

function localNeumann!(Nm::Vector, neumannFunc::Function,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)

    #∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #Nm::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        #∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dS = get_dS(element, ∂x_∂ξ, shapeFunction, ipNo)
        ϕ = get_ϕ(shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, ϕ)
        nm = neumannFunc(x; kwargs4function...)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                Nm[problemDim*(a-1)+i] += ϕ[a]*nm[i]*dS
            end
        end
    end
    return nothing
end


function localScalar!(S::Array{Float64,1}, scalarFunc::Function,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)

    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)

    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        dΩ = get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
        ϕ = get_ϕ(shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, ϕ)
        s = scalarFunc(x; kwargs4function...)
        for i ∈ 1:problemDim
            S[i] += s[i]*dΩ
        end
    end
    return nothing
end

function localScalarNeumann!(S::Array{Float64,1}, scalarFunc::Function,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)

    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ =get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        #∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dS = get_dS(element, ∂x_∂ξ, shapeFunction, ipNo)
        ϕ = get_ϕ(shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, ϕ)
        s = scalarFunc(x; kwargs4function...)
        for i ∈ 1:problemDim
            S[i] += s[i]*dS
        end
    end
    return nothing
end
