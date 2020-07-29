using RapidFEM, SparseArrays, WriteVTK

function addCell!(cells::Array{MeshCell,1}, element::LineElement)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_CURVE, element.nodeTags))
end

function addCell!(cells::Array{MeshCell,1}, element::TriElement)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_TRIANGLE, element.nodeTags))
end

function addCell!(cells::Array{MeshCell,1}, element::QuadElement)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL, element.nodeTags))
end

function addCell!(cells::Array{MeshCell,1}, element::TetElement)
    nodeTagsArray::Array{Int64} = Array{Int64,1}(undef, element.noOfElementNodes)
    if element.order ==1
        nodeTagsArray = element.nodeTags[1:4]
    elseif element.order == 2
        nodeTagsArray = [element.nodeTags[1:8]..., element.nodeTags[10], element.nodeTags[9]]
    elseif element.order == 3
        nodeTagsArray = [element.nodeTags[1:10]..., element.nodeTags[12:-1:11]..., element.nodeTags[16:-1:13]..., element.nodeTags[18], element.nodeTags[20], element.nodeTags[19], element.nodeTags[17]]
    end
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON, nodeTagsArray))
end

function addCell!(cells::Array{MeshCell,1}, element::HexElement)
    push!(cells, MeshCell(VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON, element.nodeTags))
end

function WriteToVTK(x::Vector, fieldName::String, mesh::Mesh, attributeArray::Array{Tuple{Int64, Int64},1}, problemDim::Int64)
    cells::Array{MeshCell,1} = []
    for attribute ∈ attributeArray
        for element ∈ mesh.Elements[attribute...]
            addCell!(cells, element)
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
    vtkfile = vtk_grid(fieldName, pointData, cells)
    vtkfile[fieldName] = x
    outfiles = vtk_save(vtkfile)
end

function poissonEquation()
    mesh::Mesh = RapidFEM.readMesh("../test/Tetrahedral.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 3
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,1)
    dirchAttrib::Tuple{Int64, Int64} = (2,2)
    K::SparseMatrixCSC = RapidFEM.assembleMatrix(volAttrib, FeSpace, mesh, RapidFEM.local_lagrange_K, problemDim)
    source(x) = [0.0, 0.0, 0.0]
    f::Vector = RapidFEM.assembleVector(source, volAttrib, FeSpace, mesh, RapidFEM.localSource, problemDim)
    neumann(x) = [0.0, 0.1, 0.0]
    f += RapidFEM.assembleVector(neumann, neumAttrib, FeSpace, mesh, RapidFEM.localNeumann, problemDim)
    DirichletFunction(x) = zeros(problemDim)
    RapidFEM.applyDirichletBC!(K, f, DirichletFunction, dirchAttrib, mesh, problemDim)
    x::Vector = K\f
    WriteToVTK(x, "field",mesh, [volAttrib], problemDim)
end
