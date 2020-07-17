include("gaussLine.jl")
include("gaussTri.jl")
include("gaussQuad.jl")
include("gaussTet.jl")
include("gaussHex.jl")
include("../FEM/elements.jl")

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
    return getQuadratureQuad(element.order)
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
    return getQuadratureHex(element.order)
end
