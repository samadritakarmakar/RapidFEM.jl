mutable struct VTKMeshData
    fileName::String
    mesh::Mesh
    pointData::Array{Float64, 2}
    cells::Array{MeshCell,1}
    vtkCollectionFile::WriteVTK.CollectionFile
end

function getNodeTagArray(element::AbstractElement, tagArray::Array{Int64,1})::Array{Int64}
    nodeTagsArray::Array{Int64} = Array{Int64,1}(undef, element.noOfElementNodes)
    vtkTag::Int64 = 1
    for tag ∈ tagArray
        nodeTagsArray[vtkTag] = element.nodeTags[tag]
        vtkTag +=1
    end
    return nodeTagsArray
end

function addCell!(cells::Array{MeshCell,1}, element::LineElement, tagArray::Array{Int64,1}, elementNo::Int64)
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    #push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_CURVE, element.nodeTags))
    cells[elementNo] = MeshCell(VTKCellTypes.VTK_LAGRANGE_CURVE, nodeTagsArray)
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::TriElement, tagArray::Array{Int64,1}, elementNo::Int64)
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    #push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_TRIANGLE, element.nodeTags))
    cells[elementNo] = MeshCell(VTKCellTypes.VTK_LAGRANGE_TRIANGLE, element.nodeTags, nodeTagsArray)
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::QuadElement, tagArray::Array{Int64,1}, elementNo::Int64)
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    #push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL, nodeTagsArray))
    cells[elementNo] = MeshCell(VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL, nodeTagsArray)
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::TetElement, tagArray::Array{Int64,1}, elementNo::Int64)
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    #push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON, nodeTagsArray))
    cells[elementNo] = MeshCell(VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON, nodeTagsArray)
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::HexElement, tagArray::Array{Int64,1}, elementNo::Int64)
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    #push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON, nodeTagsArray))
    cells[elementNo] = MeshCell(VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON, nodeTagsArray)
    return nothing
end

function getNodeTagDict()::Dict{Tuple{DataType, Int64, String}, Array{Int64}}
    nodeTagDict::Dict{Tuple{DataType, Int64, String}, Array{Int64}} = Dict()
    nodeTagDict[LineElement, 1, "gmsh"] = [collect(1:2)...]
    nodeTagDict[LineElement, 2, "gmsh"] = [collect(1:3)...]
    nodeTagDict[LineElement, 3, "gmsh"] = [collect(1:4)...]
    nodeTagDict[TriElement, 1, "gmsh"] = [collect(1:3)...]
    nodeTagDict[TriElement, 2, "gmsh"] = [collect(1:6)...]
    nodeTagDict[TriElement, 3, "gmsh"] = [collect(1:10)...]
    nodeTagDict[QuadElement, 1, "gmsh"] = [collect(1:4)...]
    nodeTagDict[QuadElement, 2, "gmsh"] = [collect(1:9)...]
    nodeTagDict[QuadElement, 3, "gmsh"] = [collect(1:8)..., 10, 9, 12, 11, 13, 14, 16, 15]
    nodeTagDict[TetElement, 1, "gmsh"] = [collect(1:4)...]
    nodeTagDict[TetElement, 2, "gmsh"] = [collect(1:8)..., 10, 9]
    nodeTagDict[TetElement, 3, "gmsh"] = [collect(1:10)..., 12, 11, collect(16:-1:13)..., 18, 20, 19, 17]
    nodeTagDict[HexElement, 1, "gmsh"] = [collect(1:8)...]
    nodeTagDict[HexElement, 2, "gmsh"] = [collect(1:9)..., 12, 14, 10, 17, 19, 20, 18, 11, 13, 16, 15, 23, 24, 22, 25, 21, 26, 27]
    nodeTagDict[HexElement, 3, "gmsh"] = [collect(1:10)..., 15, 16, 20, 19, 11, 12, 25, 26, 29, 30, 32, 31, 27, 28, 13, 14, 17, 18, 23, 24, 21, 22, 41, 44, 42, 43, 45, 46, 48, 47, 37, 38, 40, 39, 50, 49, 51, 52, 33, 36, 34, 35, 53, 54, 56, 55, 57, 58, 60, 59, 61, 62, 64, 63]
    return nodeTagDict
end

function createFileNameTree(filename::String, step::Int64=1)
    mkpath(filename*"/Data")
    filename *= "/Data/Data"*string(step)
    return filename
end

function InitializeVTK_Collection(fileName::String)::WriteVTK.CollectionFile
    return WriteVTK.paraview_collection(fileName)
end

function InitializeVTK(x::Vector, fileName::String, mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64},1}, problemDim::Int64)::VTKMeshData
    endElementNo::Array{Int64,1} = Array{Int64,1}(undef, length(attributeArray))
    attributeNo::Int64 = 1

    lastElement::Int64 = 0
    for attribute ∈ attributeArray
        endElementNo[attributeNo] = length(mesh.Elements[attribute...])
        attributeNo += 1
    end
    cells::Array{MeshCell,1} = Array{MeshCell,1}(undef, sum(endElementNo))
     nodeTagDict::Dict{Tuple{DataType, Int64, String}, Array{Int64}} = getNodeTagDict()
     elementNo::Int64 = 1
     attributeNo = 1
    for attribute ∈ attributeArray
        Threads.@threads for elementNo ∈ 1:endElementNo[attributeNo]
            element::AbstractElement = mesh.Elements[attribute...][elementNo]
            tagArray::Array{Int64,1} = nodeTagDict[typeof(element), element.order, mesh.meshSoftware]
            addCell!(cells, element, tagArray, lastElement+elementNo)
        end
        lastElement = endElementNo[attributeNo]
        attributeNo += 1
    end
    Nodes::Dict{Any,Any} = mesh.Nodes
    noOfNodes::Int64 = mesh.noOfNodes
    pointData::Array{Float64, 2} = Array{Float64, 2}(undef, 3, noOfNodes)
    Threads.@threads for nodeNo ∈ 1:noOfNodes
        pointData[1:3, nodeNo] = Nodes[nodeNo]'
    end
    return VTKMeshData(fileName, mesh, pointData, cells, InitializeVTK_Collection(fileName*"/"*fileName))
end

function createVTKFile(vtkMeshData::VTKMeshData, step::Int64=1)
    fileName::String = createFileNameTree(vtkMeshData.fileName, step)
    vtkfile::WriteVTK.DatasetFile = WriteVTK.vtk_grid(fileName, vtkMeshData.pointData, vtkMeshData.cells)
    return vtkfile
end

function vtkSave(vtkMeshData::VTKMeshData)
    vtkSave(vtkMeshData.vtkCollectionFile)
end

function vtkSave(vtkfile::T) where T
    WriteVTK.vtk_save(vtkfile)
    return nothing
end

function vtkDataAdd(vtkMeshData::VTKMeshData, dataTuple::T1, dataNameTuple::T2, time::Float64 = 0.0, step::Int64=1) where {T1, T2}
    @assert (length(dataTuple) == length(dataNameTuple)) "The dataTuple and dataNameTuple lengths must match."
    vtkFile::WriteVTK.DatasetFile = createVTKFile(vtkMeshData, step)
    for i ∈ 1:length(dataNameTuple)
        vtkFile[dataNameTuple[i]] = dataTuple[i]
    end
    vtkMeshData.vtkCollectionFile[time] = vtkFile
    return nothing
end
