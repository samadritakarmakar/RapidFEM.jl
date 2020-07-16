include("../Mesh/mesh.jl")

mesh = readMesh("Bar2.msh");
println("Attributes are:\n", mesh.attributes)
for attribute ∈ mesh.attributes
    println("No of Elements in attribute: ", attribute, " = ", getNoOfElements(mesh.Elements[attribute]))
end
println("Meshing Software: ", mesh.meshSoftware)
