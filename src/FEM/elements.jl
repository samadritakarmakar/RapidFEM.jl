"""This defines the an abstract type of element, to be used for functions which
are common to all kinds of elements"""
abstract type AbstractElement end

"""This defines a struct for Point Elements"""
struct PointElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

"""This defines a struct for Line Elements"""
struct LineElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

"""This defines a struct for Triangular Elements"""
struct TriElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

"""This defines a struct for Quadrilateral Elements"""
struct QuadElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

"""This defines a struct for Tetrahedral Elements"""
struct TetElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end

"""This defines a struct for Hexahedral Elements"""
struct HexElement <: AbstractElement
    label::Int64
    attributes::Array{Int64}
    nodeTags::Array{Int64}
    noOfElementNodes::Int64
    order::Int64
end
