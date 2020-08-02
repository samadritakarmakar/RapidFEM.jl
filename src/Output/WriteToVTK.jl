function getNodeTagArray(element::AbstractElement, tagArray::Array{Int64,1})::Array{Int64}
    nodeTagsArray::Array{Int64} = Array{Int64,1}(undef, element.noOfElementNodes)
    vtkTag::Int64 = 1
    for tag ∈ tagArray
        nodeTagsArray[vtkTag] = element.nodeTags[tag]
        vtkTag +=1
    end
    return nodeTagsArray
end

function addCell!(cells::Array{MeshCell,1}, element::LineElement, tagArray::Array{Int64,1})
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_CURVE, element.nodeTags))
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::TriElement, tagArray::Array{Int64,1})
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_TRIANGLE, element.nodeTags))
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::QuadElement, tagArray::Array{Int64,1})
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL, nodeTagsArray))
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::TetElement, tagArray::Array{Int64,1})
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON, nodeTagsArray))
    return nothing
end

function addCell!(cells::Array{MeshCell,1}, element::HexElement, tagArray::Array{Int64,1})
    nodeTagsArray::Array{Int64} = getNodeTagArray(element, tagArray)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON, nodeTagsArray))
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

function InitializeVTK(x::Vector, fileName::String, mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64},1}, problemDim::Int64)::WriteVTK.DatasetFile
    cells::Array{MeshCell,1} = []
     nodeTagDict::Dict{Tuple{DataType, Int64, String}, Array{Int64}} = getNodeTagDict()
    for attribute ∈ attributeArray
        for element ∈ mesh.Elements[attribute...]
            tagArray::Array{Int64,1} = nodeTagDict[typeof(element), element.order, mesh.meshSoftware]
            addCell!(cells, element, tagArray)
        end
    end
    Nodes::Dict{Any,Any} = mesh.Nodes
    pointData::Array{Float64, 2} = Array{Float64, 2}(undef, 3, mesh.noOfNodes)
    for nodeNo ∈ 1:mesh.noOfNodes
        #pointData[nodeNo, 1] = 3
        #for dim ∈ 1:problemDim
            pointData[1:3, nodeNo] = Nodes[nodeNo]'
        #end
    end
    #return cells, pointData
    vtkfile::WriteVTK.DatasetFile = WriteVTK.vtk_grid(fileName, pointData, cells)
    return vtkfile
    #=vtkfile[fieldName] = x
    outfiles = WriteVTK.vtk_save(vtkfile)=#
end

function vtkSave(vtkfile::WriteVTK.DatasetFile)
    WriteVTK.vtk_save(vtkfile)
    return nothing
end
