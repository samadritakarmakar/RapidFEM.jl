#using LinearAlgebra
"""This function creates a Dict of ShapeFunctions arrays. Since the library
uses isoparametric elements, the shape functions for a certain element are
calculated just once and used all over the domain. The returned Dict variable
keeps a list of the shape functions used in the whole domian. A new Dict variable
may be generated in the below manner:

    FeSpace = RapidFEM.createFeSpace()
"""
function createFeSpace()::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}
    feSpace = Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}()
    return feSpace
end

"""This function checks for the given element if the shape function has already been generated.
If not generated, the a new set of shape functions are generated and returned to the FeSpace variable

    shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, lagrange)
"""
function feSpace!(FeSpace::Dict{Tuple{DataType, Int64, Any}, Array{ShapeFunction}}, element::AbstractElement, mesh::Mesh, elementFunction::Function)::Array{ShapeFunction}
    #element::AbstractElement = mesh.Elements[attribute][elementNo]
    typeOfElement::DataType = typeof(element)
    if (typeOfElement, element.order, elementFunction) ∉ keys(FeSpace)
        FeSpace[typeOfElement, element.order, elementFunction] = calculateShapeFunctions(element, elementFunction, mesh.meshSoftware)
    end
    return FeSpace[typeOfElement, element.order, elementFunction]
end

"""This function return an interpolated value of the position x, a set of shape function
values for certain set of integration points.

    x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
"""
function getInterpolated_x(CoordArray::Array{Float64,2}, ϕ::Array{Float64,1})::Array{Float64,1}
    return CoordArray*ϕ
end

"""As the function name suggests, it returns the value of ∂x_∂ξ

    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)

where ipNo is the integartion point number.
"""
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

"""This function returns the appropriate function to be used with the passed
on element to calculate ∂ξ_∂x

    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
"""
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
    #return ipData.w*abs(0.5*det(∂x_∂ξ_temp))
    return ipData.w*abs(det(∂x_∂ξ_temp))
end

function get_dΩ_Tet(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    ∂x_∂ξ_temp = Array{Float64,2}(undef, size(∂x_∂ξ,1)+1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[1,:] = ones(1,size(∂x_∂ξ,2))
    ∂x_∂ξ_temp[2:end,:] = ∂x_∂ξ
    #return ipData.w*abs((1.0/6.0)*det(∂x_∂ξ_temp))
    return ipData.w*abs(det(∂x_∂ξ_temp))
end

"""This function return the appropriate function to be used to calculate
the dΩ for the given element

    dΩFunc::Function = getFunction_dΩ(element)
"""
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
    #return ipData.w*0.5*norm(cross(∂x_∂ξ[:,1],∂x_∂ξ[:,2]))
    crossProd::Array{Float64,1} = cross(∂x_∂ξ[:,2]-∂x_∂ξ[:,1],∂x_∂ξ[:,3]-∂x_∂ξ[:,1])
    return ipData.w*norm(crossProd)

end

function get_dS_Normalized(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    return ipData.w*norm(cross(∂x_∂ξ[:,1],∂x_∂ξ[:,2]))
end

function get_dL(∂x_∂ξ::Array{Float64}, ipData::IpPoint)::Float64
    return ipData.w*norm(∂x_∂ξ[:,1])
end

"""This function returns the dS for integrating over the neumann boundary

    dSFunc::Function = getFunction_dS(element)
"""
function getFunction_dS(element::TriElement)::Function
    return get_dS_Tri
end

function getFunction_dS(element::QuadElement)::Function
    return get_dS_Normalized
end

function getFunction_dS(element::LineElement)::Function
    return get_dL
end
"""This function returns the dL for integrating over a line in the boundary

    dLFunc::Function = getFunction_dL(element)
"""
function getFunction_dL(element::LineElement)::Function
    return get_dL
end
