include("lagrangeLine.jl")
include("lagrangeTri.jl")
include("lagrangeQuad.jl")
include("lagrangeTet.jl")
include("lagrangeHex.jl")
#include("../FEM/elements.jl")



function lagrange(element::LineElement, meshSoftware::String)
    if(element.order == 1 && meshSoftware == "gmsh")
        N = lagrangeLineOrder1_gmsh
    elseif(element.order == 2 && meshSoftware == "gmsh")
        N = lagrangeLineOrder2_gmsh
    elseif(element.order == 3 && meshSoftware == "gmsh")
        N = lagrangeLineOrder3_gmsh
    else
        error("Given order: "*element.order*" is not supported yet")
    end
    return N
end

function lagrange(element::TriElement, meshSoftware::String)
    if(element.order == 1 && meshSoftware == "gmsh")
        N = lagrangeTriOrder1_gmsh
    elseif(element.order == 2 && meshSoftware == "gmsh")
        N = lagrangeTriOrder2_gmsh
    elseif(element.order == 3 && meshSoftware == "gmsh")
        N = lagrangeTriOrder3_gmsh
    else
        error("Given order: "*element.order*" is not supported yet")
    end
    return N
end

function lagrange(element::QuadElement, meshSoftware::String)
    if(element.order == 1 && meshSoftware == "gmsh")
        N = lagrangeQuadOrder1_gmsh
    elseif(element.order == 2 && meshSoftware == "gmsh")
        N = lagrangeQuadOrder2_gmsh
    elseif(element.order == 3 && meshSoftware == "gmsh")
        N = lagrangeQuadOrder3_gmsh
    else
        error("Given order: "*element.order*" is not supported yet")
    end
    return N
end

function lagrange(element::TetElement, meshSoftware::String)
    if(element.order == 1 && meshSoftware == "gmsh")
        N = lagrangeTetOrder1_gmsh
    elseif(element.order == 2 && meshSoftware == "gmsh")
        N = lagrangeTetOrder2_gmsh
    elseif(element.order == 3 && meshSoftware == "gmsh")
        N = lagrangeTetOrder3_gmsh
    else
        error("Given order: "*element.order*" is not supported yet")
    end
    return N
end

function lagrange(element::HexElement, meshSoftware::String)
    if(element.order == 1 && meshSoftware == "gmsh")
        N = lagrangeHexOrder1_gmsh
    elseif(element.order == 2 && meshSoftware == "gmsh")
        N = lagrangeHexOrder2_gmsh
    elseif(element.order == 3 && meshSoftware == "gmsh")
        N = lagrangeHexOrder3_gmsh
    else
        error("Given order: "*element.order*" is not supported yet")
    end
    return N
end
