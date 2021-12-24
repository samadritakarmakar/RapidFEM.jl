using LinearAlgebra, RapidFEM
#This part of the library is based on documentation of Salome platform given in the links below
#2D : https://docs.salome-platform.org/8/gui/SMESH/aspect_ratio.html
#3D : https://docs.salome-platform.org/8/gui/SMESH/aspect_ratio_3d.html

function getConnectionData(coordArray::AbstractArray{Float64, 2}, nodeNos::Array{Array{Int64,1},1})
    connections = zeros(size(coordArray, 1), length(nodeNos))
    connectionLengths = zeros(length(nodeNos))
    connectionNo = 1
    for nodeNo ∈ nodeNos
        connections[:, connectionNo] = coordArray[:,nodeNo[2]]-coordArray[:,nodeNo[1]]
        connectionLengths[connectionNo] = norm(connections[:, connectionNo])
        connectionNo += 1
    end
    return connections, connectionLengths
end

function getSurfaceNormals(connections::Array{Float64, 2}, normalCombs::Array{Array{Int64,1},1})
    normals = zeros(size(connections, 1), length(normalCombs))
    normalMags = zeros(length(normalCombs))
    normalNo = 1
    for normalComb ∈ normalCombs
        normals[:,normalNo] = cross(connections[:,normalComb[1]], connections[:,normalComb[2]])
        normalMags[normalNo] = norm(normals[:,normalNo])
        normalNo += 1
    end
    return normals, normalMags
end

function getAspectRatioOfTriElement(coordArray::AbstractArray{Float64})
    nodeNos = [[1, 2], [1, 3], [2, 3]]
    connections, connectionLengths = getConnectionData(coordArray, nodeNos)
    A = 0.5*norm(cross(connections[:, 1], connections[:, 2]))
    connectionLengthMax = maximum(connectionLengths)
    return connectionLengthMax*sum(connectionLengths)/(4.0*sqrt(3.0)*A)
end

function getAspectRatioOfElement(mesh::Mesh, element::TriElement)
    coordArray = (getCoordArray(mesh, element))[:,1:3]
    return getAspectRatioOfTriElement(coordArray)
end
###Test with top = [0.5, 1/(2*√3), √(2/3)]; getAspectRatioOfTetElement([[0.0, 0, 0] [1, 0,0] [.5, sind(60), 0.0] top])
function getAspectRatioOfTetElement(coordArray::AbstractArray{Float64})
    nodeNos = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
    connections, connectionLengths = getConnectionData(coordArray, nodeNos)
    connectionLengthMax = maximum(connectionLengths)
    normalCombs = [[1,2],[1,3],[2,3],[4,5]] #Combination of nodeNos (connections) that create a normal each
    normals, normalMags = getSurfaceNormals(connections, normalCombs)
    α = dot(connections[:,1], normals[:,3])#normals[:,2])
    r = abs(α)/sum(normalMags)
    return connectionLengthMax/(2.0*sqrt(6.0)*r)
end

function getAspectRatioOfElement(mesh::Mesh, element::TetElement)
    coordArray = (getCoordArray(mesh, element))[:,1:4]
    return getAspectRatioOfTetElement(coordArray)
end


function getAspectRatioOfElement(mesh::Mesh, element::LineElement)
    return 1.0
end

function getAspectRatioOfElement(mesh::Mesh, element::PointElement)
    return 1.0
end

function getAspectRatioOfQuadElement(coordArray::AbstractArray{Float64})
    nodeNos = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4],[3, 4]]
    connections, connectionLengths = getConnectionData(coordArray, nodeNos)
    connectionLengthMax = maximum(connectionLengths)
    A = 0.5*norm(cross(connections[:, 1], connections[:, 2]))
    A += 0.5*norm(cross(connections[:, 4], connections[:, 6]))
    return connectionLengthMax*sum(connectionLengths)/(4*A)
end

function getAspectRatioOfElement(mesh::Mesh, element::QuadElement)
    coordArray = (getCoordArray(mesh, element))[:,1:4]
    return getAspectRatioOfQuadElement(coordArray)
end

function getAspectRatioOfElement(element::HexElement)
    coordArray = (getCoordArray(mesh, element))[:,1:8]
    error("Hexahedral Elements are not supported yet for aspect ratio.")
end


function getAspectRatios(mesh::Mesh, attribArray::Array{Tuple{Int64, Int64}, 1})
    totalElements = 0
    sort!(attribArray)
    for attrib ∈ attribArray
        totalElements +=length(mesh.Elements[attrib...])
    end
    aspectRatios = zeros(totalElements)
    elementNo = 1
    for attrib ∈ attribArray
        elements = mesh.Elements[attrib...]
        for element ∈ elements
            aspectRatios[elementNo] = getAspectRatioOfElement(mesh, element)
            elementNo += 1
        end
    end
    return aspectRatios
end