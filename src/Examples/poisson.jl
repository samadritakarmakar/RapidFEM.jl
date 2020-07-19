using RapidFEM

function local_K(element::Ab)


function poissonEquation()
    mesh = RapidFEM.readMesh("../test/Bar.msh")
    FeSpace = RapidFEM.createFeSpace()
    element = mesh.Elements[3,4][1]
    coordArray = RapidFEM.getCoordArray(mesh, element)
    shapeFunction = RapidFEM.feSpace!(FeSpace, element, mesh, RapidFEM.lagrange)
    ∂ξ_∂xFunc = RapidFEM.getFunction_∂ξ_∂x(element)
    ∂x_∂ξ = RapidFEM.get_∂x_∂ξ(coordArray, shapeFunction[1].∂ϕ_∂ξ)
    ∂ξ_dx = ∂ξ_∂xFunc(∂x_∂ξ)
    dΩFunc = getFunction_dΩ(element)
    dΩ = dΩFunc(∂x_∂ξ, shapeFunction[1].ipData)
end
