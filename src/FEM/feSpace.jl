include("../ShapeFunctions/shapeFunction.jl")
using LinearAlgebra

function createFeSpace()::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}
    feSpace = Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}()
    return feSpace
end

function feSpace!(FeSpace::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}, element::AbstractElement, mesh::Mesh, elementFunction::Function)::Array{ShapeFunction}
    #element::AbstractElement = mesh.Elements[attribute][elementNo]
    shapeFunctionAtIp = Array{ShapeFunction}(undef, 1)
    typeOfElement::DataType = typeof(element)
    if (typeOfElement, element.order) ∉ keys(FeSpace)
        FeSpace[typeOfElement, element.order] = calculateShapeFunctions(element, elementFunction, mesh.meshSoftware)
    end
    return FeSpace[typeOfElement, element.order]
end

function getInterpolated_x(CoordArray::Array{Float64,2}, ϕ::Array{Float64,1})::Array{Float64,1}
    return CoordArray*ϕ
end

function get_∂x_∂ξ(CoordArray::Array{Float64,2}, ∂ϕ_∂ξ::Array{Float64})::Array{Float64}
    return CoordArray*∂ϕ_∂ξ
end

function get_∂ξ_∂x_Normalized(∂x_∂ξ::Array{Float64})::Array{Float64}
    return inv(∂x_∂ξ)
end

function get_∂ξ_∂x_TriTet(∂x_∂ξ::Array{Float64})::Array{Float64}
    ∂x_∂ξ_temp = Array{Float64,2}(undef, size(∂x_∂ξ,1)+1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[1,:] = ones(1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[2:end,:] = ∂x_∂ξ
    ∂1_∂x__∂x_∂x::Array{Float64,2} = [zeros(size(∂x_∂ξ,1))';diagm(ones(size(∂x_∂ξ,1)))]
    return ∂x_∂ξ_temp\∂1_∂x__∂x_∂x
end

function getFunction_∂ξ_∂x(element::T)::Function where {T<:AbstractElement}
    return get_∂ξ_∂x_Normalized
end

function getFunction_∂ξ_∂x(element::TriElement)::Function
    return get_∂ξ_∂x_TriTet
end

function getFunction_∂ξ_∂x(element::TetElement)::Function
    return get_∂ξ_∂x_TriTet
end

function get_dΩ_Nomalized(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    return ipData.w*abs(det(∂x_∂ξ))
end

function get_dΩ_Tri(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64

    ∂x_∂ξ_temp = Array{Float64,2}(undef, size(∂x_∂ξ,1)+1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[1,:] = ones(1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[2:end,:] = ∂x_∂ξ
    return ipData.w*abs(0.5*det(∂x_∂ξ_temp))
end

function get_dΩ_Tet(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    ∂x_∂ξ_temp = Array{Float64,2}(undef, size(∂x_∂ξ,1)+1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[1,:] = ones(1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[2:end,:] = ∂x_∂ξ
    return ipData.w*abs((1.0/6.0)*det(∂x_∂ξ_temp))
end

function getFunction_dΩ(element::T)::Function where {T<:AbstractElement}
    return get_dΩ_Nomalized
end

function getFunction_dΩ(element::TriElement)::Function
    return get_dΩ_Tri
end

function getFunction_dΩ(element::TetElement)::Function
    return get_dΩ_Tet
end

function get_dS_Tri(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    return ipData.w*0.5*norm(cross(∂x_∂ξ[:,1],∂x_∂ξ[:,2]))
end

function get_dS_Normalized(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    return ipData.w*norm(cross(∂x_∂ξ[:,1],∂x_∂ξ[:,2]))
end

function get_dL(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    return ipData.w*norm(∂x_∂ξ[:,1])
end

function getFunction_dS(element::TriElement)::Function
    return get_dS_Tri
end

function getFunction_dS(element::QuadElement)::Function
    return get_dS_Normalized
end

function getFunction_dS(element::LineElement)::Function
    return get_dL
end

function getFunction_dL(element::LineElement)::Function
    return get_dL
end
