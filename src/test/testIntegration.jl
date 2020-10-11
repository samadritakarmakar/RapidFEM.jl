using RapidFEM
function testIntegration()
    mesh::Mesh = RapidFEM.readMesh("OneElmntMsh/TetrahedralOrder3.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 1
    volAttrib::Tuple{Int64, Int64} = (3,3)
    neumAttrib::Tuple{Int64, Int64} = (2,1)
    dirchAttrib::Tuple{Int64, Int64} = (2,2)
    activeDimensions::Array{Int64,1} = [1, 1, 1]
    #println(mesh.Elements[3,3][1].nodeTags)
    #Delete later##########################################
    element::AbstractElement = mesh.Elements[volAttrib][1]
    coordArrayTemp::Array{Float64,2} = getCoordArray(mesh, element)
    coordArray::Array{Float64,2} =  coordArrayTemp[1:1:3,:]
    #println(coordArray)
    f = open("hexList", "w")
    str = ""
    for coordNo ∈ 1:size(coordArray, 2)
        for dim ∈ 1:3
            str *= string(coordArray[dim,coordNo])*"\t"
        end
        str *= "\n"
    end
    write(f, str)
    close(f)
    scalarFunction(x) = [1.0]
    v1::Array{Float64,1} = RapidFEM.assembleScalar(scalarFunction, volAttrib, FeSpace, mesh, RapidFEM.localScalar!, problemDim, activeDimensions)
    v2::Array{Float64,1} = RapidFEM.assembleScalar(scalarFunction, neumAttrib, FeSpace, mesh, RapidFEM.localScalarNeumann!, problemDim, activeDimensions)
    return sum(v1), sum(v2)
end
