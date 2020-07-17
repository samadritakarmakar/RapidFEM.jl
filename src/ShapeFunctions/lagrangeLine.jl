function lagrangeLineOrder1_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    ξ_::Array{Float64}= [-1.0 0 0;
    1.0 0 0]
    ϕ::Array{T}= Array{T}(undef, size(ξ_, 1))
    for a ∈ 1:size(ξ_, 1)
        ξₐ::T = ξ_[a,1]
        ϕ[a] = 0.5*(1+ ξₐ*ξ)
    end
    return ϕ
end

function lagrangeLineOrder2_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    ξ_::Array{T,2} = [-1.0 0 0;
    1.0 0 0;
    0 0 0]
    ϕ::Array{T}= Array{T}(undef, size(ξ_, 1))
    for a ∈ 1:2
        ξₐ::T = ξ_[a,1]
        ϕ[a] = 0.5*(1+ ξₐ*ξ)*ξₐ*ξ
    end
    ϕ[3] = (1-ξ^2)
    return ϕ
end


function lagrangeLineOrder3_gmsh(ξ_Array::Array{T})::Array{T} where T
    ξ::T = ξ_Array[1]
    ξ_::Array{T,2} = [-1.0 0 0;
    1.0 0 0;
    -1.0/3.0 0 0;
    1.0/3.0 0 0]
    ϕ::Array{T}= Array{T}(undef, size(ξ_, 1))
    for a ∈ 1:2
        ξₐ::T = ξ_[a,1]
        ϕ[a] = -9.0/16.0*(1+ ξₐ*ξ)*(1.0/9.0-ξ^2)
    end
    for a ∈ 3:4
        ξₐ::T = ξ_[a,1]
        ϕ[a] = 27.0/16.0*(1.0/3.0+3.0*ξₐ*ξ)*(1.0-ξ^2)
    end
    return ϕ
end


function checkLagrangeLine_gmsh()
    ξ_::Array{Float64,2} = [-1.0 0 0;
    1.0 0 0;
    -1.0/3.0 0 0;
    1.0/3.0 0 0]
    for pnt ∈ 1:size(ξ_,1)
        ξ=ξ_[pnt,1]
        ϕ = lagrangeLineOrder3_gmsh(ξ)
        println(pnt)
        println("ξ= ", ξ)
        for i ∈ 1:length(ϕ)
            print(ϕ[i], " ")
        end
        println()
        println()
    end
end
