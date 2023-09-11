#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
"""This function finds the gradient of a solution u at a certain (Integration) point, given that you know the solution at the nodes of the element
using u_Nodes::AbstractArray{Float64, 1} = getSolAtElement(u, element, problemDim) and the gradient of the shape function at the point using 
∂ϕ_∂x = get_∂ϕ_∂x(element, ∂x_∂ξ, shapeFunction, ipNo)

∂u_∂x is of size (problemDim,  size(∂ϕ_∂x, 2))

    ∂u_∂x = get_∂u_∂x!(∂u_∂x, u_Nodes, ∂ϕ_∂x, problemDim)
"""
 function get_∂u_∂x!(∂u_∂x::AbstractArray{Float64, 2}, u_Nodes::AbstractArray{Float64, 1}, ∂ϕ_∂x::AbstractArray{Float64,2}, problemDim::Int64)
    fill!(∂u_∂x, 0.0)
    for a ∈ 1:size(∂ϕ_∂x, 1)
        for J ∈ 1:size(∂ϕ_∂x, 2)
            for i ∈ 1:problemDim
                ∂u_∂x[i,J] += ∂ϕ_∂x[a,J]*u_Nodes[problemDim*(a-1)+i]
            end
        end
    end
end


"""This function finds the gradient of a solution u at a certain (Integration) point, given that you know the solution at the nodes of the element
using u_Nodes::AbstractArray{Float64, 1} = getSolAtElement(u, element, problemDim) and the gradient of the shape function at the point using 
∂ϕ_∂x = get_∂ϕ_∂x(element, ∂x_∂ξ, shapeFunction, ipNo)

∂u_∂x is of size (problemDim,  size(∂ϕ_∂x, 2))

    ∂u_∂x = get_∂u_∂x(u_Nodes, ∂ϕ_∂x, problemDim)"""

function get_∂u_∂x(u_Nodes::AbstractArray{Float64, 1}, ∂ϕ_∂x::AbstractArray{Float64,2}, problemDim::Int64)
    ∂u_∂x = zeros(problemDim, size(∂ϕ_∂x, 2))
    get_∂u_∂x!(∂u_∂x, u_Nodes, ∂ϕ_∂x, problemDim)
    return ∂u_∂x
end

"""Initializes ∂u_∂x to zeros(problemDim, size(∂ϕ_∂x, 2))"""
initialize_∂u_∂x(∂ϕ_∂x::AbstractArray{Float64,2}, problemDim::Int64) = zeros(problemDim, size(∂ϕ_∂x, 2))


"""This function finds the solution u at a certain (Integration) point, given that you know the solution at the nodes of the element
using u_Nodes::AbstractArray{Float64, 1} = getSolAtElement(u, element, problemDim) and the shape function at the point using 
ϕ::AbstractArray{Float64,1} = RapidFEM.get_ϕ(shapeFunction, ipNo)

u is of size problemDim

get_u!(u::AbstractArray{Float64, 1}, u_Nodes::AbstractArray{Float64, 1}, ϕ::AbstractArray{Float64,1}, problemDim::Int64)"""

function get_u!(u::AbstractArray{Float64, 1}, u_Nodes::AbstractArray{Float64, 1}, ϕ::AbstractArray{Float64,1}, problemDim::Int64)
    fill!(u, 0.0)
    for a ∈ 1:length(ϕ)
        for i ∈ 1:problemDim
            u[i] += ϕ[a]*u_Nodes[problemDim*(a-1)+i]
        end
    end
end

"""This function finds the solution u at a certain (Integration) point, given that you know the solution at the nodes of the element
using u_Nodes::AbstractArray{Float64, 1} = getSolAtElement(u, element, problemDim) and the shape function at the point using 
ϕ::AbstractArray{Float64,1} = RapidFEM.get_ϕ(shapeFunction, ipNo)

u is of size problemDim

u  = get_u(u_Nodes::AbstractArray{Float64, 1}, ϕ::AbstractArray{Float64,1}, problemDim::Int64)"""

function get_u(u_Nodes::AbstractArray{Float64, 1}, ϕ::AbstractArray{Float64,1}, problemDim::Int64)
    u = zeros(problemDim)
    get_u!(u, u_Nodes, ϕ, problemDim)
    return u
end

"""getCurrentCoordArray returns the current coord matrix of the element given that you know the solution at the nodes of the element
using u_Nodes::AbstractArray{Float64, 1} = getSolAtElement(u, element, problemDim) and the coord matrix of the element from 

Represents the continuum concept
x = u + X

    getCurrentCoordArray!(currentCoordArray, coordArray, u_Nodes)    
"""
function getCurrentCoordArray!(currentCoordArray::AbstractArray{Float64}, coordArray::AbstractArray{Float64}, u_Nodes::AbstractArray{Float64, 1})
    problemDim = size(coordArray, 1)
    noOfNodes = size(coordArray, 2)
    for a ∈ 1:noOfNodes
        for i ∈ 1:problemDim
            currentCoordArray[i, a] = coordArray[i, a] + u_Nodes[problemDim*(a-1)+i]
        end
    end
end



"""getCurrentCoordArray returns the current coord matrix of the element given that you know the solution at the nodes of the element
using u_Nodes::AbstractArray{Float64, 1} = getSolAtElement(u, element, problemDim) and the coord matrix of the element from 

Represents the continuum concept
x = u + X

    currentCoordArray = getCurrentCoordArray(coordArray, u_Nodes)    
"""
function getCurrentCoordArray(coordArray::AbstractArray{Float64}, u_Nodes::AbstractArray{Float64, 1})
    currentCoordArray::Array{Float64} = zeros(size(coordArray)...)
    getCurrentCoordArray!(currentCoordArray, coordArray, u_Nodes)
    return currentCoordArray
end
 
"""This function finds the size of the side of the element which has the lowest inclination
to the given velocity vector

    elmntSizeAlongVel(mesh, element, velocity)
"""
function elmntSizeAlongVel(mesh::Mesh, element::AbstractElement, velocity::AbstractArray{Float64,1})
    coordArrayTrans::Array{Float64,2} = getCoordArray(mesh, element)'
    noOfRows = size(coordArrayTrans, 1)
    eye::Array{Float64, 2} = Array{Float64, 2}(I, noOfRows, noOfRows)
    permMatrix::Array{Float64, 2} = zeros(noOfRows, noOfRows)
    permMatrix[1,:] = eye[3,:]
    for i ∈ 2:noOfRows
        permMatrix[i,:] = eye[i-1,:]
    end
    coordArrayPerm::Array{Float64,2} = permMatrix*coordArrayTrans
    vecMat::Array{Float64, 2} = coordArrayTrans-coordArrayPerm
    vecNorm::Array{Float64, 1} = Array{Float64, 1}(undef, size(vecMat, 1))
    for i ∈ size(vecMat,1)
        vecNorm[i] = norm(vecMat[i,:])
    end
    cosθ::Array{Float64, 1} = vecMat*velocity/(vecNorm*norm(velocity))
    maxCosθ::Float64, n::Int64 = findmax(cosθ)
    return coordArrayTrans[n, :]'
end
