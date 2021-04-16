using RapidFEM

function testNormalShapeFunction()
    meshFiles = ["OneElmntMsh/LineOrder", "OneElmntMsh/QuadrilateralOrder", "OneElmntMsh/HexahedralOrder"]
    errorflag = false
    dim = [1, 2, 3]
    for meshFileNo ∈ 1:length(meshFiles)
        meshFile = meshFiles[meshFileNo]
        for order ∈ 1:3
            testFile = meshFile*string(order)*".msh"
            mesh = readMesh(testFile)
            element = mesh.Elements[dim[meshFileNo],3][1]
            N = lagrange(element, mesh.meshSoftware)
            println("Checking : ",N, " with Test File : ", testFile)
            coordArray = RapidFEM.getCoordArray(mesh, element)
            for nodeNo ∈ 1:size(coordArray, 2)
                ϕ = N(coordArray[:, nodeNo])
                if (abs(ϕ[nodeNo]-1.0)> 1e-12)
                    errorflag = true
                end
                for ϕ_no ∈ 1:length(ϕ)
                    if (ϕ_no !=nodeNo)
                        if (abs(ϕ[ϕ_no])> 1e-12)
                            errorflag = true
                        end
                    end 
                end
            end
            if(errorflag)
                error("Check Shape Functions in : ", N)
            end
        end
    end
    if(!errorflag)
        println("All Line, Quadrilateral and Hexahedral Element checks passed!")
    end
end