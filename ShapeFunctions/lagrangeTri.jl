function lagrangeTriOrder1_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ζ::T = ξ_Array[3]
    ϕ::Array{T} = Array{T}(undef, 3)
    #ζ = 1-(ξ+η)
    ϕ[1] = ξ
    ϕ[2] = η
    ϕ[3] = ζ
    return ϕ
end

function lagrangeTriOrder2_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ζ::T = ξ_Array[3]
    ϕ::Array{T} = Array{T}(undef, 6)
    #ζ = 1-(ξ+η)
    ϕ[1] = (2*ξ-1)*ξ
    ϕ[2] = (2*η-1)*η
    ϕ[3] = (2*ζ-1)*ζ
    ϕ[4] = 4*ξ*η
    ϕ[5] = 4*η*ζ
    ϕ[6] = 4*ζ*ξ
    return ϕ
end


function lagrangeTriOrder3_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    η::T =  ξ_Array[2]
    ζ::T = ξ_Array[3]
    ϕ::Array{T} = Array{T}(undef, 10)
    #ζ = 1-(ξ+η)
    ϕ[1] = 0.5*(3.0*ξ-1.0)*(3.0*ξ-2.0)*ξ
    ϕ[2] = 0.5*(3.0*η-1.0)*(3.0*η-2.0)*η
    ϕ[3] = 0.5*(3.0*ζ-1.0)*(3.0*ζ-2.0)*ζ
    ϕ[4] = 9.0/2.0*ξ*η*(3.0*ξ-1.0)
    ϕ[5] = 9.0/2.0*ξ*η*(3.0*η-1.0)
    ϕ[6] = 9.0/2.0*ζ*η*(3.0*η-1.0)
    ϕ[7] = 9.0/2.0*ζ*η*(3.0*ζ-1.0)
    ϕ[8] = 9.0/2.0*ξ*ζ*(3.0*ζ-1.0)
    ϕ[9] = 9.0/2.0*ξ*ζ*(3.0*ξ-1.0)
    ϕ[10] = 27*ξ*η*ζ
    return ϕ
end

function checkLagrangeTri()
    ξ_η::Array{Float64,2} = [1.0 0.0; 0.0 1.0; 0.0 0.0;
    2.0/3.0 1.0/3.0; 1.0/3.0 2.0/3.0;
    0.0 2.0/3.0; 0.0 1.0/3.0;
    1.0/3.0 0.0; 2.0/3.0 0.0;
    1.0/3.0 1.0/3.0]
    for pnt ∈ 1:size(ξ_η,1)
        ξ=ξ_η[pnt,1]
        η=ξ_η[pnt,2]
        ζ=1-(ξ+η)
        ϕ = lagrangeTriOrder3_gmsh(ξ, η, ζ)
        println(pnt)
        println("ξ= ", ξ, " η = ", η)
        for i ∈ 1:length(ϕ)
            println(ϕ[i])
        end
        println()
    end
end
