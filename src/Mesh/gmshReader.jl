#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
 
function check4seperator(readPosition::Int64, meshData::String,
    Separator::Array{Char} = [' ', '\n'])
    flag::Bool = true
    currentChar::Char = meshData[readPosition]
    for sep ∈ Separator
        if(currentChar == sep)
            flag = false
        end
    end
    return flag
end

function getNextWord(readPosition::Int64, meshData::String,
    Separator::Array{Char}= [' ', '\n'])
    flag::Bool = true
    nextWord::Array{Char} = []
    while (meshData[readPosition] == ' ')
        readPosition += 1
    end
    while(flag)
        push!(nextWord, meshData[readPosition])
        readPosition += 1
        flag = check4seperator(readPosition, meshData, Separator)
    end
    readPosition += 1
    return nextWord, readPosition
end

function getNextLine(readPosition::Int64, meshData::String)
    return getNextWord(readPosition, meshData, ['\n'])
end

function readMeshFormat(readPosition::Int64, meshData::String)
    version::Array{Char}, readPosition = getNextWord(readPosition,meshData)
    readPosition = dumpTill(readPosition, meshData, "\$EndMeshFormat")
    return parse(Float64,String(version)), readPosition
end

function dumpTill(readPosition::Int64, meshData::String, dumpStop::String)
    dump::Array{Char}, readPosition = getNextLine(readPosition,meshData)
    while (String(dump) != dumpStop)
        dump, readPosition = getNextLine(readPosition,meshData)
    end
    return readPosition
end

function dumpB4char(readPosition::Int64, meshData::String, dumpStopChar::Char)
    while meshData[readPosition] != dumpStopChar
        readPosition += 1
    end
    return readPosition
end

function checkGmshFileVersion(version::Float64)
    if (version <2 && version>3)
        error("Version: "*string(version)*" not supported.")
    end
end

function getPhysicalGroupData!(attributes::Array{Tuple{Int64, Int64},1}, AttributeName::Dict{Any,Any}, readPosition::Int64, meshData::String)
    #Reading dimension
    #AttributeName = Dict()
    noOfAttribChar::Array{Char}, readPosition = getNextWord(readPosition, meshData)
    noOfAttrib::Int64 = parse(Int64, String(noOfAttribChar))
    for i ∈ 1:noOfAttrib
    #while(meshData[readPosition+1] != '$')
        dimChar, readPosition = getNextWord(readPosition, meshData)
        dim::Int64 = parse(Int64, String(dimChar))
        AttribChar, readPosition = getNextWord(readPosition, meshData)
        Attribute::Int64 = parse(Int64, String(AttribChar))
        AttribNameChar, readPosition = getNextWord(readPosition, meshData)
        AttribName::String = String(AttribNameChar[2:end-1])
        AttributeName[dim, Attribute] = AttribName
        push!(attributes, (dim, Attribute))
    end
    readPosition = dumpTill(readPosition, meshData, "\$EndPhysicalNames")
    return noOfAttrib, readPosition
end

function getNodesData!(Nodes::Dict{Any,Any}, readPosition::Int64, meshData::String)
    NoOfNodeChar::Array{Char}, readPosition = getNextWord(readPosition, meshData)
    noOfNodes::Int64 = parse(Int64, String(NoOfNodeChar))
    for Node ∈ 1:noOfNodes
        Coord::Array{Float64} =[]
        NodeLabelChar::Array{Char}, readPosition = getNextWord(readPosition, meshData)
        NodeLabel::Int64 = parse(Int64, String(NodeLabelChar))
        for i ∈ 1:3
            CoordChar::Array{Char}, readPosition = getNextWord(readPosition, meshData)
            push!(Coord, parse(Float64, String(CoordChar)))
        end
        #println(Coord)
        Nodes[NodeLabel] = Coord
    end
    readPosition = dumpTill(readPosition, meshData, "\$EndNodes")
    return noOfNodes, readPosition
end

function getElementTypeProperty(gmshElmentType::Int64)
    elementType::Dict{Int64, Any} = Dict()
    elementType[15] = "Point"
    elementType[1] = "Line"
    elementType[8] = "Line"
    elementType[26] = "Line"
    elementType[2] = "Tri"
    elementType[9] = "Tri"
    elementType[21] = "Tri"
    elementType[3] = "Quad"
    elementType[10] = "Quad"
    elementType[36] = "Quad"
    elementType[5] = "Hex"
    elementType[12] = "Hex"
    elementType[92] = "Hex"
    elementType[4] = "Tet"
    elementType[11] = "Tet"
    elementType[29] = "Tet"

    order::Dict{Int64, Any} = Dict()
    order[15] = 0
    order[1] = 1
    order[2] = 1
    order[3] = 1
    order[4] = 1
    order[5] = 1
    order[8] = 2
    order[9] = 2
    order[10] = 2
    order[11] = 2
    order[12] = 2
    order[21] = 3
    order[26] = 3
    order[29] = 3
    order[92] = 3
    order[36] = 3
    dim::Dict{Int64, Any} = Dict()
    dim[15] = 0
    dim[1] = 1
    dim[8] = 1
    dim[26] = 1
    dim[2] = 2
    dim[9] = 2
    dim[21] = 2
    dim[3] = 2
    dim[10] = 2
    dim[36] = 2
    dim[5] = 3
    dim[12] = 3
    dim[92] = 3
    dim[4] = 3
    dim[11] = 3
    dim[29] = 3
    noOfElmntNodes::Dict{Int64, Any} = Dict()
    noOfElmntNodes[15] = 1
    noOfElmntNodes[1] = 2
    noOfElmntNodes[8] = 3
    noOfElmntNodes[26] = 4
    noOfElmntNodes[2] = 3
    noOfElmntNodes[9] = 6
    noOfElmntNodes[21] = 10
    noOfElmntNodes[3] = 4
    noOfElmntNodes[10] = 9
    noOfElmntNodes[36] = 16
    noOfElmntNodes[5] = 8
    noOfElmntNodes[12] = 27
    noOfElmntNodes[92] = 64
    noOfElmntNodes[4] = 4
    noOfElmntNodes[11] = 10
    noOfElmntNodes[29] = 20
    if gmshElmentType ∉ keys(elementType)
        error("The gmsh element type "*string(gmshElmentType)*" is not supported yet")
    end
    return dim[gmshElmentType], elementType[gmshElmentType], order[gmshElmentType],
    noOfElmntNodes[gmshElmentType]
end


function getElementData!(Elements::Dict{Any,Any}, readPosition::Int64, meshData::String)
    noOfElementsString::Array{Char}, readPosition = getNextWord(readPosition, meshData)
    noOfElements::Int64 = parse(Int64, String(noOfElementsString))
    dim0elements::Int64 = 0
    dim1elements::Int64 = 0
    dim2elements::Int64 = 0
    dim3elements::Int64 = 0
    for ElementNo ∈ 1:noOfElements
        labelString::Array{Char}, readPosition = getNextWord(readPosition, meshData)
        label::Int64 = parse(Int64, String(labelString))
        elementTypeNoString::Array{Char}, readPosition = getNextWord(readPosition, meshData)
        elementTypeNo::Int64 = parse(Int64, String(elementTypeNoString))
        dim::Int64, elementType::String,
        order::Int64, noOfElmntNodes::Int64 =
        getElementTypeProperty(elementTypeNo)
        noOfTagsString::Array{Char}, readPosition = getNextWord(readPosition, meshData)
        noOfTags::Int64 = parse(Int64, String(noOfTagsString))
        attributes::Array{Int64, 1} = []
        for tagList ∈ 1:noOfTags
            tagString::Array{Char,1}, readPosition  = getNextWord(readPosition, meshData)
            tag::Int64 = parse(Int64, String(tagString))
            push!(attributes, tag)
        end
        nodeTags::Array{Int64} = []
        for elementNo ∈ 1:noOfElmntNodes
            nodeTagString::Array{Char,1}, readPosition =  getNextWord(readPosition, meshData)
            nodeTag::Int64 = parse(Int64, String(nodeTagString))
            push!(nodeTags, nodeTag)
        end
        if elementType == "Point"
            dim0elements +=1
            #Elements[0, dim0elements] =
            if (0, attributes[1]) ∉ keys(Elements)
                Elements[0, attributes[1]] = []
            end
            #Elements[0, attributes[1]] =
            #PointElement(label, attributes, nodeTags, noOfElmntNodes, order)
            push!(Elements[0, attributes[1]],
            PointElement(label, attributes, nodeTags, noOfElmntNodes, order))
        elseif elementType == "Line"
            dim1elements +=1
            if (1, attributes[1]) ∉ keys(Elements)
                Elements[1, attributes[1]] = []
            end
            #Elements[1, dim1elements] =
            #LineElement(label, attributes, nodeTags, noOfElmntNodes, order)
            push!(Elements[1, attributes[1]],
            LineElement(label, attributes, nodeTags, noOfElmntNodes, order))
        elseif elementType == "Tri"
            dim2elements +=1
            if (2, attributes[1]) ∉ keys(Elements)
                Elements[2, attributes[1]] = []
            end
            #Elements[2, dim2elements] =
            #TriElement(label, attributes, nodeTags, noOfElmntNodes, order)
            push!(Elements[2, attributes[1]],
            TriElement(label, attributes, nodeTags, noOfElmntNodes, order))
        elseif elementType == "Quad"
            dim2elements +=1
            if (2, attributes[1]) ∉ keys(Elements)
                Elements[2, attributes[1]] = []
            end
            #Elements[2, dim2elements] =
            #QuadElement(label, attributes, nodeTags, noOfElmntNodes, order)
            push!(Elements[2, attributes[1]],
            QuadElement(label, attributes, nodeTags, noOfElmntNodes, order))
        elseif elementType == "Hex"
            dim3elements +=1
            if (3, attributes[1]) ∉ keys(Elements)
                Elements[3, attributes[1]] = []
            end
            #Elements[3, dim3elements] =
            #HexElement(label, attributes, nodeTags, noOfElmntNodes, order)
            push!(Elements[3, attributes[1]],
            HexElement(label, attributes, nodeTags, noOfElmntNodes, order))
        elseif elementType == "Tet"
            dim3elements +=1
            if (3, attributes[1]) ∉ keys(Elements)
                Elements[3, attributes[1]] = []
            end
            #Elements[3, dim3elements] =
            #TetElement(label, attributes, nodeTags, noOfElmntNodes, order)
            push!(Elements[3, attributes[1]],
            TetElement(label, attributes, nodeTags, noOfElmntNodes, order))
        end
    end
    readPosition = dumpTill(readPosition, meshData, "\$EndElements")
    return noOfElements, readPosition
end




function readGmshFile!(gmshFileName::String, attributes::Array{Tuple{Int64, Int64},1},
    AttributeName::Dict{Any,Any}, Nodes::Dict{Any,Any}, Elements::Dict{Any, Any})
    io::IOStream = open(gmshFileName)
    meshData::String = read(io, String)
    dataLength::Int64 = length(meshData)
    readPosition::Int64 = 1
    noOfAttrib::Int64 = 0
    #attributes::Array{Tuple{Int64, Int64},1} = []
    #AttributeName::Dict{Any,Any} = Dict()
    noOfNodes::Int64 = 0
    #Nodes::Dict{Any,Any} = Dict()
    noOfElements::Int64 =0
    #Elements::Dict{Any, Any} = Dict()
    while (readPosition < dataLength)
        readPosition = dumpB4char(readPosition, meshData, '$')
        word, readPosition = getNextWord(readPosition,meshData)
        readItems = ["\$MeshFormat" "\$PhysicalNames" "\$Nodes" "\$Elements"]
        currentWord::String = String(word)
        #Reading Mesh Format
        if currentWord == readItems[1]
            version::Float64, readPosition = readMeshFormat(readPosition,meshData)
            checkGmshFileVersion(version)
        elseif currentWord == readItems[2]
            #noOfAttrib, AttributeName, readPosition = getPhysicalGroupData(readPosition, meshData)
            noOfAttrib, readPosition = getPhysicalGroupData!(attributes, AttributeName, readPosition, meshData)
        elseif currentWord == readItems[3]
            noOfNodes, readPosition = getNodesData!(Nodes, readPosition, meshData)
        elseif currentWord == readItems[4]
            noOfElements, readPosition = getElementData!(Elements, readPosition, meshData)
        end
    end
    #println(Nodes)
    #println()
    #println(Elements)
    return noOfAttrib, noOfNodes, noOfElements
end
