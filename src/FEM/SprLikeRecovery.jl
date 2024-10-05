using LinearAlgebra, RapidFEM
function getPolyTri(coord::AbstractVector{Float64}, order::Int64)
    x, y = coord[1], coord[2]
    if order == 0
        return [1.0]
    elseif order == 1
        return [1.0, x, y]
    elseif order == 2
        return [1.0, x, y, x^2, x*y, y^2]
    elseif order == 3
        return [1.0, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3]
    else
        error("Order not implemented")
    end
end

function getPolyTri(coords::AbstractMatrix{Float64}, order::Int64)
    P = zeros(Float64, size(coords, 2), (order+1)*(order+2)÷2)
    i = 1
    for coord in eachcol(coords)
        P[i, :] = getPolyTri(coord, order)
        i += 1
    end
    return P
end

getPoly(::TriElement, coords::Union{AbstractVector{Float64}, AbstractMatrix{Float64}}, order::Int64) = getPolyTri(coords, order)


function getPolyQuad(coord::AbstractVector{Float64}, order::Int64)
    x, y = coord[1], coord[2]
    poly = Array{Float64, 1}(undef, (order+1)^2)
    term = 1
    for i in 0:order
        for j in 0:order
            # Construct the term x^i * y^j
            poly[term] = x^i * y^j 
            term += 1
        end
    end
    return poly
end

function getPolyQuad(coords::AbstractMatrix{Float64}, order::Int64)
    P = zeros(Float64, size(coords, 2), (order+1)^2)
    i = 1
    for coord in eachcol(coords)
        P[i, :] = getPolyQuad(coord, order)
        i += 1
    end
    return P
end

getPoly(::QuadElement, coords::Union{AbstractVector{Float64}, AbstractMatrix{Float64}}, order::Int64) = getPolyQuad(coords, order)

function getPolyLine(coord::AbstractVector{Float64}, order::Int64)
    x = coord[1]
    poly = Array{Float64, 1}(undef, (order+1))
    term = 1
    for i in 0:order
        # Construct the term x^i * y^j
        poly[term] = x^i
        term += 1
    end
    return poly
end

function getPolyLine(coords::AbstractMatrix{Float64}, order::Int64)
    P = zeros(Float64, size(coords, 2), order+1)
    i = 1
    for coord in eachcol(coords)
        P[i, :] = getPolyLine(coord, order)
        i += 1
    end
    return P
end

getPoly(::LineElement, coords::Union{AbstractVector{Float64}, AbstractMatrix{Float64}}, order::Int64) = getPolyLine(coords, order)

function getPolyTet(coord::AbstractVector{Float64}, order::Int64)
    x, y, z = coord[1], coord[2], coord[3]
    if order == 0
        return [1.0]
    elseif order == 1
        return [1.0, x, y, z]
    elseif order == 2
        return [1.0, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2]
    elseif order == 3
        return [1.0, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2, x^3, x^2*y, x^2*z, x*y^2, x*y*z, x*z^2, y^3, y^2*z, y*z^2, z^3]
    else
        error("Order not implemented")
    end
end

function getPolyTet(coords::AbstractMatrix{Float64}, order::Int64)
    P = zeros(Float64, size(coords, 2), (order+1)*(order+2)*(order+3)÷6)
    i = 1
    for coord in eachcol(coords)
        P[i, :] = getPolyTet(coord, order)
        i += 1
    end
    return P
end

getPoly(::TetElement, coords::Union{AbstractVector{Float64}, AbstractMatrix{Float64}}, order::Int64) = getPolyTet(coords, order)

function getPolyHex(coord::AbstractVector{Float64}, order::Int64)
    x, y, z = coord[1], coord[2], coord[3]
    poly = Array{Float64, 1}(undef, (order+1)^3)
    term = 1
    for i in 0:order
        for j in 0:order
            for k in 0:order
                # Construct the term x^i * y^j * z^k
                poly[term] = x^i * y^j * z^k
                term += 1
            end
        end
    end
    return poly
end

function getPolyHex(coords::AbstractMatrix{Float64}, order::Int64)
    P = zeros(Float64, size(coords, 2), (order+1)^3)
    i = 1
    for coord in eachcol(coords)
        P[i, :] = getPolyHex(coord, order)
        i += 1
    end
    return P
end

getPoly(::HexElement, coords::Union{AbstractVector{Float64}, AbstractMatrix{Float64}}, order::Int64) = getPolyHex(coords, order)

function getTotalIpPoints(attribElementNos::Vector{Tuple{Tuple{Int64, Int64}, Int64}}, mesh::Mesh, FeSpace::Dict, reduction::Int64)
    totalIpPoints = 0
    for attribElementNo ∈ attribElementNos
        element = mesh.Elements[attribElementNo[1]][attribElementNo[2]]
        shapeFunction = feSpace!(FeSpace, element, mesh, reduction = reduction)
        totalIpPoints += length(shapeFunction)
    end
    return totalIpPoints
end

function getSprPolysAroundNode(centerNodeNo::Int64, FeSpace::Dict, mesh::Mesh, meshExtra::MeshExtra, usedDims::StepRange; 
    reduction::Int64 = 0)
    
    # get elements around centerNodeNo
    attribElementNos = meshExtra.nodeToElementMap[centerNodeNo]
    totalIpPoints = getTotalIpPoints(attribElementNos, mesh, FeSpace, reduction)
    #the type of polynomial and order to be used is as per the first element
    firstElement = mesh.Elements[attribElementNos[1][1]][attribElementNos[1][2]]
    nodeCoord = mesh.Nodes[centerNodeNo][usedDims]
    #println("totalIpPoints: ", totalIpPoints, " length(p): ", length(p))
    #=if length(p) >= totalIpPoints
        #add more neaby elements to attribElementNos
        newAttribElementNos = deepcopy(attribElementNos)
        for attribElementNo ∈ attribElementNos
            element = mesh.Elements[attribElementNo[1]][attribElementNo[2]]
            nodeTags = element.nodeTags
            for nodeTag ∈ nodeTags
                append!(newAttribElementNos, meshExtra.nodeToElementMap[nodeTag])
            end
        end
        unique!(sort!(newAttribElementNos))
        attribElementNos = newAttribElementNos
        totalIpPoints = getTotalIpPoints(attribElementNos, mesh, FeSpace, reduction)
        @assert length(p) < totalIpPoints "Error: Not enough Integration points to perform SPR"
    end=#
    sort!(attribElementNos)
    ipCoords = zeros(Float64, length(usedDims), totalIpPoints)
    totalIpNo = 1
    for attribElementNo ∈ attribElementNos
        element = mesh.Elements[attribElementNo[1]][attribElementNo[2]]
        shapeFunction = feSpace!(FeSpace, element, mesh, reduction = reduction)
        coordArray = getCoordArray(mesh, element)[usedDims, :]
        for ipNo in 1:length(shapeFunction)
            ϕ = get_ϕ(shapeFunction, ipNo)
            ipCoords[:, totalIpNo] = getInterpolated_x(coordArray, ϕ)
            totalIpNo += 1
        end
    end
    orderUsed = firstElement.order
    p = getPoly(firstElement, nodeCoord, orderUsed)
    while length(p) > totalIpPoints && orderUsed > 0
        orderUsed -= 1
        p = getPoly(firstElement, nodeCoord, orderUsed)
    end
    #println("orderUsed: ", orderUsed)
    @assert orderUsed >= 0 "Error: Selected order is less than 0"
    P = getPoly(firstElement, ipCoords, orderUsed)
    return p, P, attribElementNos, totalIpPoints
end

function getSprSampleMatrix(ipDataDict::Dict{Tuple{Int64, Int64}, T}, problemDim::Int64, attribElementNos::Vector{Tuple{Tuple{Int64, Int64}, Int64}}, totalIpPoints::Int64) where T
    
    sampleMatrix = Array{Float64, 2}(undef, totalIpPoints, problemDim)

    ipNo = 1
    for attribElementNo ∈ attribElementNos
        elementNo = attribElementNo[2]
        while haskey(ipDataDict, (elementNo, ipNo))
            sampleMatrix[ipNo, :] = vec(ipDataDict[elementNo, ipNo])
            ipNo += 1
        end
        ipNo = 1
    end
    return sampleMatrix
end

function SprLikeRecovery(ipDataDict::Dict, FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}}, 
    mesh::Mesh,  attrib::Tuple{Int64, Int64}, problemDim::Int64, activeDimensions::Array{Int64,1}=[1, 1, 1]; 
    reduction::Int64 = 0)

    meshExtra = MeshExtra(mesh, [attrib])
    dimRange = RapidFEM.createDimRange()
    usedDims = dimRange[activeDimensions]
    recoveredData = Array{Float64, 1}(undef, length(mesh.Nodes)*problemDim)
    for node ∈ enumerate(mesh.Nodes)
        nodeNo, nodeCoord = node
        #println("nodeNo: ", nodeNo)
        p, P, attribElementNos, totalIpPoints = getSprPolysAroundNode(nodeNo, FeSpace, mesh, meshExtra, usedDims, reduction = reduction)
        sampleMatrix = getSprSampleMatrix(ipDataDict, problemDim, attribElementNos, totalIpPoints)
        #println("size of sampleMatrix: ", size(sampleMatrix))
        #display(sampleMatrix)
        #println("size of P: ", size(P))
        recoveredData[problemDim*(nodeNo-1)+1:problemDim*nodeNo] = p'*(P \ sampleMatrix) 
    end
    return recoveredData
end

function SprLikeRecovery(postProcessFunction::func, sol::AbstractVector{Float64},
    parameters::param,  FeSpace::Dict{Tuple{DataType, Int64, Any, Int64}, Array{ShapeFunction}},
    mesh::Mesh,  attrib::Tuple{Int64, Int64}, problemDim::Int64,
    activeDimensions::Array{Int64,1}=[1, 1, 1]) where {func, param}

    dimRange = RapidFEM.createDimRange()
    usedDims = dimRange[activeDimensions]
    meshExtra = MeshExtra(mesh, [attrib])
    ipDataDict = Dict{Tuple{Int64, Int64}, AbstractArray}()
    elementNo = 1
    for element ∈ mesh.Elements[attrib]
        shapeFunction = feSpace!(FeSpace, element, mesh)
        coordArray = getCoordArray(mesh, element)
        postPrcssElIpData = postProcessFunction(parameters, sol, problemDim, element, elementNo, shapeFunction, coordArray)
        for ipNo ∈ 1:length(shapeFunction)
            if postPrcssElIpData isa Array{Float64, 2}
                ipDataDict[elementNo, ipNo] = vec(postPrcssElIpData[ipNo, :])
            elseif postPrcssElIpData isa Array{Float64, 3}
                ipDataDict[elementNo, ipNo] = vec(postPrcssElIpData[ipNo, :, :])
            elseif postPrcssElIpData isa Array{Float64, 5}
                ipDataDict[elementNo, ipNo] = vec(postPrcssElIpData[ipNo, :, :, :, :])
            else
                ipDataDict[elementNo, ipNo] = postPrcssElIpData[ipNo]
            end
        end
        elementNo += 1
    end
    recoveredData = SprLikeRecovery(ipDataDict, FeSpace, mesh, attrib, problemDim, activeDimensions)
    return recoveredData
end
            
    




        








    