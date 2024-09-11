function addGmshPhysicalNames!(returnChar::Array{Char, 1}, mesh::RapidFEM.Mesh)
    attribNames = mesh.AttributeName
    append!(returnChar, collect("\$PhysicalNames\n"))
    append!(returnChar, collect(string(length(keys(attribNames)))*"\n"))
    attribString = ""
    for attrib ∈ sort!(collect(keys(attribNames)))
        attribString *= string(attrib[1])*""" """*string(attrib[2])*""" """*"\""*string(attribNames[attrib])*"\""*"\n"
    end
    append!(returnChar, attribString)
    append!(returnChar, collect("\$EndPhysicalNames\n"))
end
function addGmshNodes!(returnChar::Array{Char, 1}, mesh::RapidFEM.Mesh)
    nodesKeyArray = sort(collect(keys(mesh.Nodes)))
    noOfNodes = length(nodesKeyArray)
    append!(returnChar, collect("\$Nodes\n$noOfNodes\n"))
    for nodeKey ∈ nodesKeyArray
        append!(returnChar, collect("$nodeKey"))
        for coordDim ∈ 1:3
            append!(returnChar, collect(" $((mesh.Nodes[nodeKey])[coordDim])"))
        end
        append!(returnChar, '\n')
    end
    append!(returnChar, "\$EndNodes\n")
end

function getGmshElementTypeNo(element::RapidFEM.PointElement)
    return 15
end

function getGmshElementTypeNo(element::RapidFEM.LineElement)
    if element.order == 1
        return 1
    elseif element.order == 2
        return 8
    elseif element.order == 3
        return 26
    end
end

function getGmshElementTypeNo(element::RapidFEM.TriElement)
    if element.order == 1
        return 2
    elseif element.order == 2
        return 9
    elseif element.order == 3
        return 21
    end
end

function getGmshElementTypeNo(element::RapidFEM.QuadElement)
    if element.order == 1
        return 3
    elseif element.order == 2
        return 10
    elseif element.order == 3
        return 36
    end
end

function getGmshElementTypeNo(element::RapidFEM.TetElement)
    if element.order == 1
        return 4
    elseif element.order == 2
        return 11
    elseif element.order == 3
        return 29
    end
end

function getGmshElementTypeNo(element::RapidFEM.HexElement)
    if element.order == 1
        return 5
    elseif element.order == 2
        return 12
    elseif element.order == 3
        return 92
    end
end

function addGmshElements!(returnChar::Array{Char, 1}, mesh::RapidFEM.Mesh)
    append!(returnChar, collect("\$Elements\n$(getTotalElements(mesh))\n"))
    elAttribs = sort(collect(keys(mesh.Elements)))
    elNo = 1
    for attrib ∈ elAttribs
        for element ∈ mesh.Elements[attrib...]
            append!(returnChar, collect("$elNo $(getGmshElementTypeNo(element))"))
            append!(returnChar, collect(" $(length(attrib))"))
            #for attribPart ∈ attrib
            #    append!(returnChar, collect(" $attribPart"))
            #end
            append!(returnChar, collect(" $(attrib[2]) $(attrib[2])"))
            for nodeTag ∈ element.nodeTags
                append!(returnChar, collect(" $nodeTag"))
            end
            elNo += 1
            append!(returnChar, '\n')
        end
    end
    append!(returnChar, collect("\$EndElements\n"))
end

function getGmshFileString(mesh::RapidFEM.Mesh)
    returnChar = collect("\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n")
    addGmshPhysicalNames!(returnChar, mesh)
    addGmshNodes!(returnChar, mesh)
    addGmshElements!(returnChar, mesh)
    return String(returnChar)
end

function writeMesh(fileName::String, mesh::RapidFEM.Mesh)
    if mesh.meshSoftware == "gmsh"
        io = open(fileName, "w")
        gmshString = getGmshFileString(mesh::RapidFEM.Mesh)
        write(io, gmshString)
        close(io)
    end
end
