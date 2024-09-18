"""Adds new boundary to the mesh given a condition on the boundarySelectCondition is satisfied, 
which takes in an argument of Vector{Float64}"""
function addBoundaryElements!(mesh::Mesh, 
    meshExtra::MeshExtra, newAttribute::Tuple, newAttributeName::String, boundarySelectCondition::Function, usedDims::StepRange)

    satisfyingNodes = Set{Int64}()
    #get nodes that satisfy boundary 
    for boundaryNode ∈ meshExtra.boundaryNodes
        nodeCoord = mesh.Nodes[boundaryNode][usedDims]
        if boundarySelectCondition(nodeCoord)
            push!(satisfyingNodes, boundaryNode)
        end
    end
    #find faces which have all nodes satisfying boundary condition
    boundaryFaces = findall(x->length(x) == 1, meshExtra.allFaces)
    #println("boundaryFaces = $boundaryFaces")
    for boundaryFace ∈ boundaryFaces
        if length(setdiff(boundaryFace, satisfyingNodes))  == 0 #complete face on boundary
            attribElNo =  meshExtra.allFaces[boundaryFace][1]
            element = mesh.Elements[attribElNo[1]][attribElNo[2]]
            faceCombinations = RapidFEM.getFaceCombinations(typeof(element), element.order)
            for faceCombination ∈ faceCombinations #one of the face combinations is a genuine boundary element
                
                faceNodes = element.nodeTags[faceCombination]
                if length(setdiff(faceNodes, boundaryFace)) == 0
                    if newAttribute ∉ keys(mesh.Elements)
                        mesh.Elements[newAttribute] = Array{AbstractElement, 1}()
                        mesh.noOfAttrib += 1
                        mesh.AttributeName[newAttribute] = newAttributeName
                        push!(mesh.attributes, newAttribute)
                    end
                    mesh.noOfElements += 1
                    #println("faceNodes $faceNodes")
                    push!(mesh.Elements[newAttribute], createNewElement(newAttribute, 1, faceNodes))
                end
            end
        end
    end
end

function addBoundaryElements!(mesh::Mesh, 
    meshExtra::MeshExtra, newAttribute::Tuple, 
    newAttributeName::String, boundarySelectCondition::Function, activeDims::Vector{Int64})

    dimRange = createDimRange()
    return addBoundaryElements!(mesh, meshExtra, newAttribute, newAttributeName, boundarySelectCondition, dimRange[activeDims])
end

"""Adds new boundaries to the mesh given the condition on the boundarySelectConditions are satisfied. Each function of boundarySelectCondition
 takes in an argument of Vector{Float64}"""
function addBoundaryElements!(mesh::Mesh, boundaryDimension::Int64, newAttributes::Vector{Tuple}, 
    newAttributeNames::Vector{String}, boundarySelectConditions::Vector{Function}, activeDims::Vector{Int64})
    dimRange = createDimRange()
    meshExtra = MeshExtra(mesh, boundaryDimension+1)
    attribNo = 1
    for newAttribute ∈ newAttributes
        addBoundaryElements!(mesh, meshExtra, newAttribute, newAttributeNames[attribNo], 
        boundarySelectConditions[attribNo], dimRange[activeDims])
    end
end
