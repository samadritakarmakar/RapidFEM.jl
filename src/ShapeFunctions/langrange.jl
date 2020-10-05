include("lagrangeLine.jl")
include("lagrangeTri.jl")
include("lagrangeQuad.jl")
include("lagrangeTet.jl")
include("lagrangeHex.jl")
#include("../FEM/elements.jl")


"""This function returns the kind of function to be used to generate the
shape function values for the given element. It also differentiates as per the
order of the element and the meshing software used. It is assumed that the meshing
software used is the same for meshing all over the domain.
"""
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
