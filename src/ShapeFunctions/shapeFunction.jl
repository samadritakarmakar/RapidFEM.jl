#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
 
include("langrange.jl")
include("../Quadrature/quadrature.jl")
include("../Mesh/mesh.jl")

using ForwardDiff

"""This is a struct that contains the weights and the integration points.
You would find this being used as a part of array of the ShapeFunction type
array in which the data can typically be accessed in the below way:

    shapeFunction[ipNo].ipData.w
    shapeFunction[ipNo].ipData.ip

where ipNo is the integartion point number.
"""
struct IpPoint
    w::Float64
    ip::Array{Float64}
end

"""This is a struct that contains the shape functions, its gradients and its hessian.
You would find this being used as a part of array of the ShapeFunction type
array in which the data can typically be accessed in the below way:

    shapeFunction[ipNo].ϕ
    shapeFunction[ipNo].∂ϕ_∂ξ
    shapeFunction[ipNo].∂²ϕ_∂ξ²

where ipNo is the integartion point number.
"""
struct ShapeFunction
    ϕ::Array{Float64,1}
    ∂ϕ_∂ξ::Array{Float64,2}
    ∂²ϕ_∂ξ²::Array{Float64,3}
    ipData::IpPoint
end

function vector_hessian(f::Function, x::Array{Float64,1}, length_f::Int64)::Array{Float64,3}
    n::Int64 = length(x)
    out::Array{Float64} = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(f, x), x)
    return reshape(out, length_f, n, n)
end

"""This function is responsible for calculation of the shape function array, its gradients and
the hessians at all the given integration points for that order of the element.

    calculateShapeFunctions(element, elementFunction, meshSoftware)
"""
function calculateShapeFunctions(element::T, elementFunction::Function, meshSoftware::String)::Array{ShapeFunction} where {T<:AbstractElement}
    w::Array{Float64}, ip::Array{Float64} = getQuadrature(element)
    shapeFunctionAtIp::Array{ShapeFunction} = []
    N::Function = elementFunction(element, meshSoftware)
    for ipNo ∈ 1:length(w)
        ϕ = N(ip[ipNo,:])
        x = ip[ipNo,:]
        ipData::IpPoint = IpPoint(w[ipNo], x)
        ∂ϕ_∂ξ::Array{Float64,2} =convert(Array{Float64,2}, ForwardDiff.jacobian(N, x))
        #∂²ϕ_∂ξ² = ForwardDiff.hessian(N, x)
        ∂²ϕ_∂ξ²::Array{Float64,3} =convert(Array{Float64,3}, vector_hessian(N, x, length(ϕ)))
        push!(shapeFunctionAtIp, ShapeFunction(ϕ, ∂ϕ_∂ξ, ∂²ϕ_∂ξ², ipData))
    end
    return shapeFunctionAtIp
end
