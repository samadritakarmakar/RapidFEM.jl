abstract type AbstractElement end

struct PointElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

struct LineElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

struct TriElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

struct QuadElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

struct TetElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

struct HexElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end
