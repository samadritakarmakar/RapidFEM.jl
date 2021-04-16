using RapidFEM
function testIntegration()
    meshFiles = ["OneElmntMsh/LineOrder", "OneElmntMsh/TriangleOrder","OneElmntMsh/QuadrilateralOrder", "OneElmntMsh/TetrahedralOrder",
    "OneElmntMsh/HexahedralOrder"]
    errorflag = false
    dim = [1, 2, 2, 3, 3]
    testVol = [2.0, 0.5, 4, 1.0/6, 8]
    activeDimensionsArray::Array{Array{Int64,1}} = [[1, 0, 0],[1, 1, 0],[1, 1, 1]]
    FeSpace = RapidFEM.createFeSpace()
    for meshFileNo ∈ 1:length(meshFiles)
        meshFile = meshFiles[meshFileNo]
        volAttrib::Tuple{Int64, Int64} = (dim[meshFileNo],3)
        neumAttrib::Tuple{Int64, Int64} = (dim[meshFileNo],1)
        dirchAttrib::Tuple{Int64, Int64} = (dim[meshFileNo],2)
        activeDimensions = activeDimensionsArray[dim[meshFileNo]]
        for order ∈ 1:3
            testFile = meshFile*string(order)*".msh"
            mesh = readMesh(testFile)
            problemDim::Int64 = 1
            #println(mesh.Elements[3,3][1].nodeTags)
            #Delete later##########################################
            #=element::AbstractElement = mesh.Elements[volAttrib][1]
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
            close(f)=#
            scalarFunction(x) = [1.0]
            println("Checking for integration against file :", testFile)
            v1::Array{Float64,1} = RapidFEM.assembleScalar(scalarFunction, volAttrib, FeSpace, mesh, RapidFEM.localScalar!, problemDim, activeDimensions)
            #v2::Array{Float64,1} = RapidFEM.assembleScalar(scalarFunction, neumAttrib, FeSpace, mesh, RapidFEM.localScalarNeumann!, problemDim, activeDimensions)
            if !(abs(sum(v1)-testVol[meshFileNo])<1e-14)
                error("Check Integration for element and order type related to file : ", testFile)
                errorflag = true
            end
        end
    end
    if !errorflag
        println("Integration test passed!")
    end
end
