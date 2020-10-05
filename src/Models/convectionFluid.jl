"""Vector type variant of the convection equation. Here the velocity function, λ, can be made dependent on the position, x

    local_v_λ_∇u_Vector(velocityFunction, problemDim, element, shapeFunction, coordArray)
"""
function local_v_λ_∇u_Vector(velocityFunction::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})::Array{Float64,2}
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    K::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        ϕ::Array{Float64} = shapeFunction[ipNo].ϕ
        x::Array{Float64, 1} = getInterpolated_x(coordArray, ϕ)
        λ::Array{Float64, 1} = parameterFunction(x)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        for b ∈ 1:noOfNodes
            for a ∈ 1:noOfNodes
                for j ∈ 1:size(∂ϕ_∂x,2)
                    for i ∈ 1:problemDim
                        K[problemDim*(a-1)+j,problemDim*(b-1)+i] += ϕ[a]*λ[i]*∂ϕ_∂x[b,j]*dΩ
                        #println("K ", problemDim*(a-1)+i," ,", problemDim*(b-1)+i, " = ", λ[i]*∂ϕ_∂x[a,j]*∂ϕ_∂x[b,j]*dΩ)
                    end
                end
            end
        end
    end
    return K
end

"""Scalar type variant of the convection equation. Here the velocity function, λ, can be made dependent on the position, x

    local_v_λ_∇u_Scalar(velocityFunction, problemDim, element, shapeFunction, coordArray)
"""
function local_v_λ_∇u_Scalar(velocityFunction::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})::Array{Float64,2}
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    K::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        ϕ::Array{Float64} = shapeFunction[ipNo].ϕ
        x::Array{Float64, 1} = getInterpolated_x(coordArray, ϕ)
        λ::Array{Float64, 1} = parameterFunction(x)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        for b ∈ 1:noOfNodes
            for a ∈ 1:noOfNodes
                for j ∈ 1:size(∂ϕ_∂x,2)
                    K[a,b] += ϕ[a]*λ[j]*∂ϕ_∂x[b,j]*dΩ
                    #println("K ", problemDim*(a-1)+i," ,", problemDim*(b-1)+i, " = ", λ[i]*∂ϕ_∂x[a,j]*∂ϕ_∂x[b,j]*dΩ)
                end
            end
        end
    end
    return K
end
