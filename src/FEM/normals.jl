#using RapidFEM
function getElementSurfaceNormal(element::LineElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64}, 
    internalPoint::Union{Vector{Float64}, Nothing} = nothing)

    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    normals = zeros(size(coordArray, 1), noOfIpPoints)
    for ipNo ∈ 1:noOfIpPoints
        ϕ = get_ϕ(shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, ϕ)
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        normals[:, ipNo] = normalize([∂x_∂ξ[2], -∂x_∂ξ[1]])
        if !isnothing(internalPoint)
            if dot(normals[:, ipNo], x-internalPoint) < 0.0
                normals[:, ipNo] *= -1.0 
            end
        end
    end
    return normals
end

function getElementSurfaceNormal(element::TriElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64}, 
    internalPoint::Union{Vector{Float64}, Nothing} = nothing)

    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    normals = zeros(size(coordArray, 1), noOfIpPoints)
    for ipNo ∈ 1:noOfIpPoints
        ϕ = get_ϕ(shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, ϕ)
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        normals[:, ipNo] = normalize(cross(∂x_∂ξ[:,3]-∂x_∂ξ[:,1], ∂x_∂ξ[:,1]-∂x_∂ξ[:,2]))
        if !isnothing(internalPoint)
            if dot(normals[:, ipNo], x-internalPoint) < 0.0
                normals[:, ipNo] *= -1.0 
            end
        end
    end
    return normals
end


function getElementSurfaceNormal(element::QuadElement, elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64},
    internalPoint::Union{Vector{Float64}, Nothing} = nothing)

    noOfIpPoints = getNoOfElementIpPoints(shapeFunction)
    noOfNodes = getNoOfElementNodes(shapeFunction)
    normals = zeros(size(coordArray, 1), noOfIpPoints)
    for ipNo ∈ 1:noOfIpPoints
        ϕ = get_ϕ(shapeFunction, ipNo)
        x = getInterpolated_x(coordArray, ϕ)
        ∂x_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction, ipNo)
        normals[:, ipNo] = normalize(cross(∂x_∂ξ[:,1],-∂x_∂ξ[:,2]))
        if !isnothing(internalPoint)
            if dot(normals[:, ipNo], x-internalPoint) < 0.0
                normals[:, ipNo] *= -1.0 
            end
        end
    end
    return normals
end

function getInternalPointOfElement(element::Union{TriElement, QuadElement, TetElement, HexElement}, 
    elementNo::Int64, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64})
    ipNo = 1
    ϕ = get_ϕ(shapeFunction, ipNo)
    x = getInterpolated_x(coordArray, ϕ)
    return x
end

function getInternalPoints(surfAttrib::Tuple{Int64, Int64}, mesh::Mesh, meshExtra::MeshExtra, FeSpace::Dict, activeDims::Vector{Int64}, dimRange::Dict{Vector{Int64}, StepRange{Int64, Int64}},
    u::Union{Nothing, Array{Float64}}= nothing, problemDim::Int64 = 0; elementFunction::Function=lagrange, quadrature::Function = gauss)

    #dimRange = RapidFEM.createDimRange()
    usedDims = dimRange[activeDims]
    internalPoints = Dict{Int64, Vector{Float64}}()
    for (elementNo, element) ∈ enumerate(mesh.Elements[surfAttrib])
        associatedVolElements = meshExtra.allFaces[sort(element.nodeTags)]
        @assert length(associatedVolElements) == 1 "Element with node tags: $(element.nodeTags) is not a boundary element"
        volElement_attrib_elNo = associatedVolElements[1]
        volElement = mesh.Elements[volElement_attrib_elNo[1]][volElement_attrib_elNo[2]]
        if  volElement isa TetElement || volElement isa TriElement
            if volElement.order == 2
                reduction = 1
            elseif volElement.order == 3
                reduction = 2
            else
                reduction = 0
            end
        elseif volElement isa HexElement || volElement isa QuadElement
            if volElement.order == 2
                reduction = 2
            elseif volElement.order == 3
                reduction = 3
            else
                reduction = 1
            end
        end
        shapeFunction = feSpace!(FeSpace, volElement, mesh, reduction = reduction, elementFunction = elementFunction, quadrature = quadrature)
        coordArray = getCoordArray(mesh, volElement)[usedDims, :]
        if !isnothing(u)
            @assert problemDim != 0 "problemDim must be set if u is not nothing"
            u_Nodes = getSolAtElement(u, volElement, problemDim)
            coordArray = getCurrentCoordArray(coordArray, u_Nodes)
        end
        internalPoints[elementNo] = getInternalPointOfElement(volElement, volElement_attrib_elNo[2], shapeFunction, coordArray)
    end
    return internalPoints
end

"""Outputs the surface normals of the elements in surfAttrib. The surface normals are checked for direction and are flipped if they point inwards.
This checking is done by using the corresponding volume element used to declare meshExtra. u is is displacement field, to be used for current configuration problems.
problemDim is the dimension of the problem. reduction is for reduced integration. elementFunction is the element function to be used. quadrature is the quadrature function to be used."""
function getSurfaceNormals(surfAttrib::Tuple{Int64, Int64}, mesh::Mesh, meshExtra::MeshExtra, FeSpace::Dict, activeDims::Vector{Int64},
    u::Union{Nothing, Array{Float64}}= nothing, problemDim::Int64 = 0; reduction = 0, elementFunction::Function=lagrange, quadrature::Function = gauss)

    dimRange = RapidFEM.createDimRange()
    usedDims = dimRange[activeDims]
    internalPoints = getInternalPoints(surfAttrib, mesh, meshExtra, FeSpace, activeDims, dimRange, u, problemDim, elementFunction = elementFunction, quadrature = quadrature)
    normals = Dict{Int64, Array{Float64}}()
    for (elementNo, element) ∈ enumerate(mesh.Elements[surfAttrib])
        shapeFunction = feSpace!(FeSpace, element, mesh, reduction = reduction, elementFunction = elementFunction, quadrature = quadrature)
        coordArray = getCoordArray(mesh, element)[usedDims, :]
        if !isnothing(u)
            @assert problemDim != 0 "problemDim must be set if u is not nothing"
            u_Nodes = getSolAtElement(u, element, problemDim, activeDims)
            coordArray = getCurrentCoordArray(coordArray, u_Nodes)
        end
        normals[elementNo] = getElementSurfaceNormal(element, elementNo, shapeFunction, coordArray, internalPoints[elementNo])
    end
    return normals
end





