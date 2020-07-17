include("gmshReader.jl")

mutable struct Mesh
    noOfAttrib::Int64
    attributes::Array{Tuple{Int64, Int64},1}
    AttributeName::Dict{Any, Any}
    noOfNodes::Int64
    Nodes::Dict{Any,Any}
    noOfElements::Int64
    Elements::Dict{Any,Any}
    meshSoftware::String
end

function readMesh(meshFileName::String)
    attributes::Array{Tuple{Int64, Int64},1} = []
    AttributeName::Dict{Any,Any} = Dict()
    Nodes::Dict{Any,Any} = Dict()
    Elements::Dict{Any, Any} = Dict()
    noOfAttrib::Int64 = 0
    noOfNodes::Int64 = 0
    noOfElements::Int64 = 0
    meshSoftware::String = ""
    posExnt = findlast('.', meshFileName)
    if meshFileName[posExnt:end] == ".msh"
        noOfAttrib, noOfNodes, noOfElements = readGmshFile!(meshFileName, attributes,
            AttributeName, Nodes, Elements)
            meshSoftware = "gmsh"
    else
        error("This mesh type is not yet supported.")
    end
    mesh::Mesh = Mesh(noOfAttrib, attributes, AttributeName, noOfNodes, Nodes, noOfElements, Elements, meshSoftware)
    return mesh
end

function getNoOfElements(ElementsArray::Array{Any,1})
    return length(ElementsArray)
end

function getCoordArray(mesh::Mesh, attribute::Tuple{Int64, Int64}, elementNo::Int64)::Array{Float64,2} 
    nodeTags::Array{Int64} = mesh.Elements[attribute][elementNo].nodeTags
    noOfElementNodes::Int64 = mesh.Elements[attribute][elementNo].noOfElementNodes
    CoordArray::Array{Float64,2} = Array{Float64}(undef, 3, noOfElementNodes)
    for elmntNodeNum::Int64 ∈ 1:noOfElementNodes
        CoordArray[:,elmntNodeNum] = mesh.Nodes[nodeTags[elmntNodeNum]]
    end
    return CoordArray
end
