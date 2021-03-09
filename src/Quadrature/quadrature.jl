#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#

include("gaussTri.jl")
include("gaussTet.jl")
include("gaussHex.jl")

function getQuadrature(element::LineElement)
    return getQuadratureLine(element.order)
end

function getQuadrature(element::TriElement)
    w::Array{Float64}, ip_base::Array{Float64} = getQuadratureTri(element.order)
    ip::Array{Float64} = Array{Float64}(undef, size(ip_base,1), 3)
    ip[:, 1:2] = ip_base
    for i ∈ 1:size(ip,1)
        ip[i,3] = 1.0 - sum(ip_base[i,:])
    end
    return w, ip
end

function getQuadrature(element::QuadElement)
    #order::Int64 = 2*element.order
    order::Int64 = Int64(ceil((element.order +1)/2))
    return getQuadratureQuad(order)
end

function getQuadrature(element::TetElement)
    w::Array{Float64}, ip_base::Array{Float64} = getQuadratureTet(element.order)
    ip::Array{Float64} = Array{Float64}(undef, size(ip_base,1), 4)
    ip[:, 1:3] = ip_base
    for i ∈ 1:size(ip,1)
        ip[i,4] = 1.0 - sum(ip_base[i,:])
    end
    return w, ip
end

function getQuadrature(element::HexElement)
    #order::Int64 = 3*element.order
    order::Int64 = Int64(ceil((element.order +1)/2)+2)
    return getQuadratureHex(order)
end
