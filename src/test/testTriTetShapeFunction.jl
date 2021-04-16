include("cart2baryTri.jl")
include("cart2baryTet.jl")

function testTriTetShapeFunction()
    meshFiles = ["OneElmntMsh/TriangleOrder", "OneElmntMsh/TetrahedralOrder"]
    errorflag = false
    dim = [2, 3]
    for meshFileNo ∈ 1:length(meshFiles)
        meshFile = meshFiles[meshFileNo]
        for order ∈ 1:3
            testFile = meshFile*string(order)*".msh"
            mesh = readMesh(testFile)
            element = mesh.Elements[dim[meshFileNo],3][1]
            N = lagrange(element, mesh.meshSoftware)
            println("Checking : ",N, " with Test File : ", testFile)
            coordArray = getCoordArray(mesh, element)
            if meshFileNo ==  1
                baryCoords = getTriBaryCentricCoordSet(Array(coordArray'))'
            elseif meshFileNo == 2
                baryCoords = getTetBaryCentricCoordSet(Array(coordArray'))'
            end
            for nodeNo ∈ 1:size(baryCoords, 2)
                ϕ = N(baryCoords[:, nodeNo])
                if (abs(ϕ[nodeNo]-1.0)> 1e-13)
                    errorflag = true
                end
                for ϕ_no ∈ 1:length(ϕ)
                    if (ϕ_no !=nodeNo)
                        if (abs(ϕ[ϕ_no])> 1e-13)
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
        println("All Triangle and Tetetrahedral Element ShapeFunction checks passed!")
    end
end