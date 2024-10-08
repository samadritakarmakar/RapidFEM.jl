#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
 
function lagrangeQuadOrder1_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ξ_η::Array{T} = [-1. -1.; 1. -1.; 1. 1.; -1. 1.]
    ϕ::Array{T}= Array{T}(undef, size(ξ_η, 1))
    for a ∈ 1:size(ξ_η, 1)
        ξ_a::T = ξ_η[a,1]
        η_a::T = ξ_η[a,2]
        ϕ[a] = (0.25*(1+ ξ_a*ξ)*(1+η_a*η))
    end
    return ϕ
end

function lagrangeQuadOrder2_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ξ_η::Array{T} = [-1. -1.; 1. -1.; 1. 1.; -1. 1.; 0 -1; 1. 0.; 0. 1.; -1. 0.; 0. 0.]
    ϕ::Array{T}= Array{T}(undef, size(ξ_η, 1))
    ϕ[1] = ((ξ - 1.0) * (ξ - 0.0) * (η - 1.0) * (η - 0.0) / 4.0)
    ϕ[2] = ((ξ + 1.0) * (ξ - 0.0) * (η - 1.0) * (η - 0.0) / 4.0)
    ϕ[3] = ((ξ + 1.0) * (ξ - 0.0) * (η + 1.0) * (η - 0.0) / 4.0)
    ϕ[4] = ((ξ - 1.0) * (ξ - 0.0) * (η + 1.0) * (η - 0.0) / 4.0)
    ϕ[5] = ((ξ + 1.0) * (ξ - 1.0) * (η - 1.0) * (η - 0.0) / -2.0)
    ϕ[6] = ((ξ + 1.0) * (ξ - 0.0) * (η + 1.0) * (η - 1.0) / -2.0)
    ϕ[7] = ((ξ - 1.0) * (ξ + 1.0) * (η + 1.0) * (η - 0.0) / -2.0)
    ϕ[8] = ((ξ - 1.0) * (ξ - 0.0) * (η + 1.0) * (η - 1.0) / -2.0)
    ϕ[9] = ((ξ - 1.0) * (ξ - -1.0) * (η + 1.0) * (η - 1.0) / 1.0)
    return ϕ
end

function lagrangeQuadOrder3_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ξ_η::Array{T} = [-1. -1.; 1. -1.; 1. 1.; -1. 1.;
                          -1.0/3. -1.; 1.0/3. -1.; 1. -1.0/3.; 1. 1.0/3.;
                            1.0/3. 1.; -1.0/3. 1.;-1. 1.0/3.; -1. -1.0/3.;
                            -1.0/3. -1.0/3.; 1.0/3. -1.0/3.; 1.0/3. 1.0/3.; -1.0/3. 1.0/3.]
    ϕ::Array{T}= Array{T}(undef, size(ξ_η, 1))
    for a ∈ 1:4
        ξ_a::T = ξ_η[a,1]
        η_a::T = ξ_η[a,2]
        ϕ[a] = (81.0/256.0*(1.0+ξ_a*ξ)*(1.0+η_a*η)*(1.0/9.0-ξ^2)*(1.0/9.0-η^2))
    end
    for a ∈ [5 6 9 10]
        ξ_a::T = ξ_η[a,1]
        η_a::T = ξ_η[a,2]
        #quad_order_3_at_5_6_9_10!(ϕ[a], ξ_a, η_a)\
        ϕ[a] = (243.0/256.0*(1.0-ξ^2)*(η^2-1.0/9.)*(1.0/3.0+3.0*ξ_a*ξ)*(1.0+η_a*η))
    end
    for a ∈ [7 8 11 12]
        ξ_a::T = ξ_η[a,1]
        η_a::T = ξ_η[a,2]
        #quad_order_3_at_7_8_11_12!(ϕ[a], ξ_a, η_a)
        ϕ[a] = (243.0/256.0*(1.0-η^2)*(ξ^2-1.0/9.)*(1.0/3.0+3.0*η_a*η)*(1.0+ξ_a*ξ))
    end
    for a ∈ 13:16
        ξ_a::T = ξ_η[a,1]
        η_a::T = ξ_η[a,2]
        ϕ[a] = (729.0/256.0*(1.0-ξ^2)*(1.0-η^2)*(1.0/3.0+3.0*ξ_a*ξ)*(1.0/3.0+3.0*η_a*η))
    end
    return ϕ
end

function checkLagrangeQuad()
    ξ_η = [-1. -1.; 1. -1.; 1. 1.; -1. 1.; 0 -1; 1. 0.; 0. 1.; -1. 0.; 0. 0.]

    for pnt ∈ 1:size(ξ_η,1)
        ξ=ξ_η[pnt,1]
        η=ξ_η[pnt,2]
        ϕ = lagrangeQuadOrder2_gmsh([ξ η])
        println(pnt)
        println("ξ= ", ξ, " η = ", η)
        for i ∈ 1:length(ϕ)
            println(ϕ[i])
        end
        println()
    end
end
