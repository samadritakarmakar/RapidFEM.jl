using RapidFEM, LinearAlgebra

defaultExtrudeFunction(x::Vector) = x

"""distributionFactor dictates how the layers are distrbuted. 
If distributionFactor = 1.0, the layers are distributed equally.
If distributionFactor > 1.0, the layers are desnser at the bottom.
If distributionFactor < 1.0, the layers are denser at the top."""
function layerDistribution(extrudeLength::Number, noOfLayers::Int, distributionFactor::Float64= 1.0)
    m = extrudeLength/(noOfLayers^distributionFactor)
    zLayers = m*collect(1:noOfLayers).^distributionFactor
end

function tapereExtrusionFunction(x::Vector, startRadius::Number, endRadius::Number, taperLength::Number) 
    x1, x2, x3 = x
    R = norm([x1, x2])
    θ = R<1e-8 ? 0.0 : atan(x2/x1)
    Rnew = R/startRadius*((1.0-x3/taperLength)*startRadius + endRadius*x3/taperLength)
    xNew = Rnew*cos(θ)
    yNew = Rnew*sin(θ)
    return [xNew, yNew, x3]
end

"""Returns a new linear hexahedral mesh with the surface extruded to a volume.
WORKS ONLY WITH LINEAR QUADRILATERAL SURFACES!
extrudeSurfaceAttrib is the attribute of the surface to be extruded.
noOfLayers is the number of layers to extrude. 
newExtrudedAttrib is the attribute of the new extruded volume. attribName is the name of the new attribute.
By default, the extrusion is parallel to the z-axis, but a custom extrusion function can be provided. For example, to extrude a surface with tapering, one can
use the tapereExtrusionFunction and define it as

    extrusionFunction(x) = tapereExtrusionFunction(x, startRadius, endRadius, taperLength).

    The distributionFactor dictates how the layers are distrbuted. 
    If distributionFactor = 1.0, the layers are distributed equally.
    If distributionFactor > 1.0, the layers are denser at the bottom.
    If distributionFactor < 1.0, the layers are denser at the top.
    
    extrudeMeshSurfaceToVolume(mesh, (2,0), 15.0, 12, (3,3), 
    "vol"; extrusionFunction=extrusionFunction, distributionFactor=0.85)
"""
function extrudeMeshSurfaceToVolume(mesh::Mesh, extrudeSurfaceAttrib::Tuple, extrudeLength::Number, noOfLayers::Int, 
    newExtrudedAttrib::Tuple, attribName::String; extrusionFunction::Function=defaultExtrudeFunction, distributionFactor::Float64=1.0)

    #zLayers = LinRange(0.0, extrudeLength, noOfLayers+1)[2:end]
    zLayers = layerDistribution(extrudeLength, noOfLayers, distributionFactor)
    meshCopy = deepcopy(mesh)

    meshCopy.Elements = Dict{Tuple{Int64, Int64}, Array{AbstractElement, 1}}()
    meshCopy.Elements[newExtrudedAttrib] = Array{AbstractElement, 1}(undef, 0)
    meshCopy.attributes = [newExtrudedAttrib]
    meshCopy.AttributeName = Dict{Tuple{Int64, Int64}, String}(newExtrudedAttrib=>attribName)

    lastNodeLayer = keys(mesh.Nodes)
    NodesCopy = deepcopy(mesh.Nodes)
    originalNodeLength = length(lastNodeLayer)
    lastElementSurface = Matrix{Int64}(undef, 4, length(mesh.Elements[extrudeSurfaceAttrib]))
    nextElementSurface = Matrix{Int64}(undef, 4, length(mesh.Elements[extrudeSurfaceAttrib]))
    for (i, element) in enumerate(mesh.Elements[extrudeSurfaceAttrib])
        @assert length(element.nodeTags) == 4 "Only linear quadrilateral elements are supported."
        lastElementSurface[:, i] = element.nodeTags
    end
    elementNo = 1
    for (layerNo, z) in enumerate(zLayers)
        #println("Layer: ", layerNo, " z: ", z)
        mapNodes = Dict{Int,Int}()
        nextNodeLayer = Int64[]
        for (i, nodeTag) in enumerate(lastNodeLayer)
            #Every node in the last layer is mapped to the next layer. An old node is mapped to a new node.
            push!(nextNodeLayer, i + layerNo*originalNodeLength)
            mapNodes[nodeTag] = nextNodeLayer[end]
            #println("Node: ", nodeTag, " nextNode: ", nextNodeLayer[end])
            if layerNo == 1
                NodesCopy[nextNodeLayer[end]] = [meshCopy.Nodes[nodeTag][1], meshCopy.Nodes[nodeTag][2]]
            else
                NodesCopy[nextNodeLayer[end]] = [NodesCopy[nodeTag][1], NodesCopy[nodeTag][2]]
            end
            
            meshCopy.Nodes[nextNodeLayer[end]] = extrusionFunction([NodesCopy[nodeTag][1], NodesCopy[nodeTag][2], z])
        end
        lastNodeLayer = deepcopy(nextNodeLayer)
        for i in 1:size(lastElementSurface, 2)
            newNodeTags =[lastElementSurface[1, i], lastElementSurface[2, i], lastElementSurface[3, i], lastElementSurface[4, i], 
            mapNodes[lastElementSurface[1, i]], mapNodes[lastElementSurface[2, i]], mapNodes[lastElementSurface[3, i]], mapNodes[lastElementSurface[4, i]]]
            nextElementSurface[:, i] = [mapNodes[lastElementSurface[1, i]], mapNodes[lastElementSurface[2, i]], mapNodes[lastElementSurface[3, i]], mapNodes[lastElementSurface[4, i]]]
            push!(meshCopy.Elements[newExtrudedAttrib], createNewElement(newExtrudedAttrib, elementNo, newNodeTags))
            elementNo += 1
       end
         lastElementSurface = deepcopy(nextElementSurface)
    end
    meshCopy.noOfNodes = length(meshCopy.Nodes)
    meshCopy.noOfElements = getTotalElements(meshCopy)
    #writeMesh("extrudedMesh.msh", meshCopy)
    return meshCopy
end
