function local_∇v_λ_∇u!(K::Array{Float64,2}, parameters::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #K[:,:] = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
        λ::Array{Float64, 1} = parameters(x)
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
    return nothing
end

function local_v_ρ_u!(M::Array{Float64,2}, parameters::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #M::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        ϕ::Array{Float64} = shapeFunction[ipNo].ϕ
        x::Array{Float64, 1} = getInterpolated_x(coordArray, ϕ)
        ρ::Array{Float64, 1} = parameters(x)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
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

function localBoundary_v_ρ_u!(M::Array{Float64,2}, parameters::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #M::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        ϕ::Array{Float64} = shapeFunction[ipNo].ϕ
        x::Array{Float64, 1} = getInterpolated_x(coordArray, ϕ)
        ρ::Array{Float64, 1} = parameters(x)
        dS::Float64 = dSFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
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


function localSource!(S::Vector, sourceFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #S::Vector = zeros(noOfNodes*problemDim)
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
    return nothing
end

function localNeumann!(Nm::Vector, neumannFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    #∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #Nm::Vector = zeros(noOfNodes*problemDim)
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
    return nothing
end


function localScalar!(S::Array{Float64,1}, scalarFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #S::Array{Float64,1} = zeros(problemDim)
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
    return nothing
end

function localScalarNeumann!(S::Array{Float64,1}, scalarFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    #∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #S::Array{Float64,1} = zeros(problemDim)
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
    return nothing
end
