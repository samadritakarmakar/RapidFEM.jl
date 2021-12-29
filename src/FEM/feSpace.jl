#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
 
"""This function creates a Dict of ShapeFunctions arrays. Since the library
uses isoparametric elements, the shape functions for any type of element is
calculated just once and used all over the domain. The returned Dict variable
keeps a list of the shape functions, it's gradients and hessians used in the whole domian. A new Dict variable
may be generated in the below manner:

    FeSpace = RapidFEM.createFeSpace()
"""
function createFeSpace()::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}
    feSpace = Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}()
    return feSpace
end

"""This function checks for the given element if the shape function has already been generated.
If not generated, the a new set of shape functions are generated and returned to the FeSpace variable

    shapeFunction::Array{ShapeFunction,1} = feSpace!(FeSpaceThreaded[currentThread], element, mesh, reduction = 0, elementFunction=lagrange)
"""
function feSpace!(FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, element::AbstractElement, mesh::Mesh; reduction::Int64 = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)::Array{ShapeFunction}
    #element::AbstractElement = mesh.Elements[attribute][elementNo]
    typeOfElement::DataType = typeof(element)
    if (typeOfElement, element.order, elementFunction) ∉ keys(FeSpace)
        FeSpace[typeOfElement, element.order, elementFunction, reduction] = calculateShapeFunctions(element, mesh.meshSoftware; elementFunction = elementFunction, reduction = reduction, quadrature = quadrature)
    end
    return FeSpace[typeOfElement, element.order, elementFunction, reduction]
end

"""This function return an interpolated value of the position x in the element. This is the value of x at the integration point.
A set of shape functions at a certain integration point is used to calculate it.

    x::Array{Float64, 1} = getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
"""
function getInterpolated_x(CoordArray::Array{Float64,2}, ϕ::Array{Float64,1})::Array{Float64,1}
    return CoordArray*ϕ
end

function getInterpolated_x(coordArray::Array{Float64,2}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)::Array{Float64,1}
    return getInterpolated_x(coordArray, shapeFunction[ipNo].ϕ)
end

"""As the function name suggests, it returns the value of ∂x_∂ξ

    ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)

where ipNo is the integartion point number.
"""
function get_∂x_∂ξ(CoordArray::Array{Float64,2}, ∂ϕ_∂ξ::Array{Float64})::Array{Float64}
    return CoordArray*∂ϕ_∂ξ
end

function get_∂x_∂ξ(coordArray::Array{Float64,2}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)::Array{Float64}
    return get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
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

function get_∂ξ_∂x(element::Union{LineElement, QuadElement, HexElement}, ∂x_∂ξ::Array{Float64})::Array{Float64}
    return get_∂ξ_∂x_Normalized(∂x_∂ξ)
end

function get_∂ξ_∂x(element::Union{TriElement, TetElement}, ∂x_∂ξ::Array{Float64})::Array{Float64}
    return get_∂ξ_∂x_TriTet(∂x_∂ξ)
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


"""Using multiple dispatch, dΩ can now be directly calculated using get_dΩ
    
    get_dΩ(element, ∂x_∂ξ, shapeFunction, ipNo)
"""
function get_dΩ(element::Union{LineElement, QuadElement, HexElement}, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)
    return get_dΩ_Nomalized(∂x_∂ξ, shapeFunction[ipNo].ipData)
end

function get_dΩ(element::TriElement, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)
    return get_dΩ_Tri(∂x_∂ξ, shapeFunction[ipNo].ipData)
end

function get_dΩ(element::TetElement, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)
    return get_dΩ_Tet(∂x_∂ξ, shapeFunction[ipNo].ipData)
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

function getFunction_dS(element::AbstractElement)::Function
    foo_surface(∂x_∂ξ::Array{Float64}, ipData::IpPoint) = 
    error("Surface Integral for this element does not exist!")
    return foo_surface
end
"""This function returns the dL for integrating over a line in the boundary

    dLFunc::Function = getFunction_dL(element)
"""
function getFunction_dL(element::LineElement)::Function
    return get_dL
end

"""Using multiple dispatch, dS can now be directly calculated using get_dS
    
    get_dS(element, ∂x_∂ξ, shapeFunction, ipNo)
"""
function get_dS(element::QuadElement, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)
    return get_dS_Normalized(∂x_∂ξ, shapeFunction[ipNo].ipData)
end

function get_dS(element::TriElement, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)
    return get_dS_Tri(∂x_∂ξ, shapeFunction[ipNo].ipData)
end

function get_dS(element::LineElement, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)
    return get_dL(∂x_∂ξ, shapeFunction[ipNo].ipData)
end

function get_dL(element::LineElement, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64)
    get_dL(∂x_∂ξ, shapeFunction[ipNo].ipData)
end

"""This function returns the ∂ϕ_∂x based on the ∂x_∂ξ returned from the get_∂x_∂ξ function"""
function get_∂ϕ_∂x(element::T, ∂x_∂ξ::Array{Float64}, shapeFunction::Array{ShapeFunction}, ipNo::Int64) where{T<:AbstractElement}
    return shapeFunction[ipNo].∂ϕ_∂ξ*get_∂ξ_∂x(element, ∂x_∂ξ)
end

"""This returns the shape Function ϕ at a certain Integration point

    get_ϕ(shapeFunction, ipNo)
"""
get_ϕ(shapeFunction::Array{ShapeFunction}, ipNo::Int64) = shapeFunction[ipNo].ϕ

getNoOfElementIpPoints(shapeFunction::Array{ShapeFunction}) = length(shapeFunction)


getNoOfElementNodes(shapeFunction::Array{ShapeFunction}) = size(shapeFunction[1].∂ϕ_∂ξ,1) 