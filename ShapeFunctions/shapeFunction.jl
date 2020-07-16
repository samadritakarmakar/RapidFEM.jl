include("langrange.jl")
include("../Quadrature/quadrature.jl")
include("../Mesh/mesh.jl")

using ForwardDiff

struct ShapeFunction
    ϕ::Array{Float64,1}
    ∂ϕ_∂ξ::Array{Float64,2}
    ∂²ϕ_∂ξ²::Array{Float64,3}
end

function vector_hessian(f::Function, x::Array{Float64,1}, length_f::Int64)::Array{Float64,3}
    n::Int64 = length(x)
    out::Array{Float64} = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(f, x), x)
    return reshape(out, length_f, n, n)
end

function calculateShapeFunctions(element::T, elementFunction::Function, meshSoftware::String)::Array{ShapeFunction} where {T<:AbstractElement}
    w::Array{Float64}, ip::Array{Float64} = getQuadrature(element)
    shapeFunctionAtIp::Array{ShapeFunction} = []
    N::Function = elementFunction(element, meshSoftware)
    for ipNo ∈ 1:length(w)
        ϕ = N(ip[ipNo,:])
        x = ip[ipNo,:]
        ∂ϕ_∂ξ::Array{Float64,2} =convert(Array{Float64,2}, ForwardDiff.jacobian(N, x))
        #∂²ϕ_∂ξ² = ForwardDiff.hessian(N, x)
        ∂²ϕ_∂ξ²::Array{Float64,3} =convert(Array{Float64,3}, vector_hessian(N, x, length(ϕ)))
        push!(shapeFunctionAtIp, ShapeFunction(ϕ, ∂ϕ_∂ξ, ∂²ϕ_∂ξ²))
    end
    return shapeFunctionAtIp
end
