#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
 
function lagrangeTetOrder1_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ζ::T = ξ_Array[3]
    γ::T = ξ_Array[4]
    ϕ::Array{T}= Array{T}(undef, 4)
    #γ = 1-(ξ+η+ζ)
    ϕ[1] = ξ
    ϕ[2] = η
    ϕ[3] = ζ
    ϕ[4] = γ
    return ϕ
end

function lagrangeTetOrder2_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ζ::T = ξ_Array[3]
    γ::T = ξ_Array[4]
    #γ = 1-(ξ+η+ζ)
    ϕ::Array{T}= Array{T}(undef, 10)
    array_ξηζγ::Array{T} = [ξ η ζ γ]
    for a ∈ 1:4
        ϕ[a] = (2.0*array_ξηζγ[a]-1.0)*array_ξηζγ[a]
    end
    pos::Array{Int64,2} = [1 2;
    2 3;
    1 3;
    1 4;
    3 4;
    2 4]
    for a ∈ 5:10
        ϕ[a] = 4*array_ξηζγ[pos[a-4,1]]*array_ξηζγ[pos[a-4,2]]
    end
    return ϕ
end

function lagrangeTetOrder3_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ζ::T = ξ_Array[3]
    γ::T = ξ_Array[4]
    #γ = 1-(ξ+η+ζ)
    ϕ::Array{T}= Array{T}(undef, 20)
    array_ξηζγ::Array{T} = [ξ η ζ γ]
    for a ∈ 1:4
        ϕ[a] = 0.5*(3.0*array_ξηζγ[a]-1.0)*(3.0*array_ξηζγ[a]-2.0)*array_ξηζγ[a]
    end
    for a ∈ 5:16
        pos2_3 = [1 2 2 3 3 1 4 1 4 3 4 2]
        pos1_3 = [2 1 3 2 1 3 1 4 3 4 2 4]
        ϕ[a] = 9/2*array_ξηζγ[pos2_3[a-4]]*array_ξηζγ[pos1_3[a-4]]*(3.0*array_ξηζγ[pos2_3[a-4]]-1.0)
    end
        ϕ[17] = 27*ξ*η*ζ
        ϕ[18] = 27*ξ*η*γ
        ϕ[19] = 27*ξ*ζ*γ
        ϕ[20] = 27*η*ζ*γ
    return ϕ
end

function checkLagrangeTet()
    #=ξ_η_ζ = [1.0  0.0  0.0  0.0;
 0.0  1.0  0.0  0.0;
 0.0  0.0  1.0  0.0;
 0.0  0.0  0.0  1.0;
 0.5  0.5  0.0  0.0;
 0.0  0.5  0.5  0.0;
 0.5  0.0  0.5  0.0;
 0.5  0.0  0.0  0.5;
 0.0  0.0  0.5  0.5;
 0.0  0.5  0.0  0.5]=#
 ξ_η_ζ = [1.0 0.0 0.0 0.0;
 0.0 1.0 0.0 0.0;
 0.0 0.0 1.0 0.0;
 0.0 0.0 0.0 1.0;
 2/3 1/3 0.0 0.0;
 1/3 2/3 0.0 0.0;
 0.0 2/3 1/3 0.0;
 0.0 1/3 2/3 0.0;
 1/3 0.0 2/3 0.0;
 2/3 0.0 1/3 0.0;
 1/3 0.0 0.0 2/3;
 2/3 0.0 0.0 1/3;
 0.0 0.0 1/3 2/3;
 0.0 0.0 2/3 1/3;
 0.0 1/3 0.0 2/3;
 0.0 2/3 0.0 1/3;
 1/3 1/3 1/3 0.0;
 1/3 1/3 0.0 1/3;
 1/3 0.0 1/3 1/3;
 0.0 1/3 1/3 1/3]
    for pnt ∈ 1:size(ξ_η_ζ,1)
        ξ=ξ_η_ζ[pnt,1]
        η=ξ_η_ζ[pnt,2]
        ζ=ξ_η_ζ[pnt,3]
        γ = 1-(ξ+η+ζ)
        ϕ = lagrangeTetOrder3_gmsh([ξ η ζ γ])
        println(pnt)
        println("ξ= ", ξ, " η = ", η, " ζ = ", ζ)
        for i ∈ 1:length(ϕ)
            println(ϕ[i])
        end
        println()
    end
end
