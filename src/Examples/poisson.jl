using RapidFEM, FEMSparse, SparseArrays, LinearAlgebra

function local_lagrange_K(problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = RapidFEM.getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    K::Array{Float64,2} = zeros(noOfNodes*problemDim, noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = RapidFEM.get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        for b ∈ 1:noOfNodes
            for a ∈ 1:noOfNodes
                for j ∈ 1:size(∂ϕ_∂x,2)
                    for i ∈ 1:problemDim
                        K[problemDim*(a-1)+i,problemDim*(b-1)+i] += ∂ϕ_∂x[a,j]*∂ϕ_∂x[b,j]*dΩ
                    end
                end
            end
        end
    end
    return K
end

function localSource(sourceFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = RapidFEM.getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    S::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = RapidFEM.get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ::Float64 = dΩFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        x::Array{Float64,1} = RapidFEM.getInterpolated_x(coordArray, ϕ)
        s::Vector = sourceFunc(x)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                S[problemDim*(a-1)+i] = ϕ[a]*s[i]*dΩ
            end
        end
    end
    return S
end

function localNeumann(neumannFunc::Function, problemDim::Int64, element::AbstractElement, shapeFunction::Array{ShapeFunction}, coordArray::Array{Float64,2})
    ∂ξ_∂xFunc::Function = RapidFEM.getFunction_∂ξ_∂x(element)
    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    Nm::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂x_∂ξ::Array{Float64,2} = RapidFEM.get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dS::Float64 = dSFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        x::Array{Float64,1} = RapidFEM.getInterpolated_x(coordArray, ϕ)
        nm::Vector = neumannFunc(x)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                Nm[problemDim*(a-1)+i] = ϕ[a]*nm[i]*dS
            end
        end
    end
    return Nm
end

function assembleMatrix(attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}, mesh::Mesh, localMatrixFunc::Function, problemDim::Int64)
    K_COO::SparseMatrixCOO = SparseMatrixCOO()
    #K_local::Array{Float64,2} = Array{Float64,2}(undef, 0, 0)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArray::Array{Float64,2} = RapidFEM.getCoordArray(mesh, element)
        shapeFunction::Array{ShapeFunction,1} = RapidFEM.feSpace!(FeSpace, element, mesh, RapidFEM.lagrange)
        K_local::Array{Float64,2} = localMatrixFunc(problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = RapidFEM.getVectorNodes(element, problemDim)
        FEMSparse.assemble_local_matrix!(K_COO, vNodes, vNodes, K_local)
    end
    return SparseArrays.SparseMatrixCSC(K_COO)
end

function assembleVector(problemfunction::Function, attribute::Tuple{Int64, Int64}, FeSpace::Dict{Tuple{DataType, Int64}, Array{ShapeFunction}}, mesh::Mesh, localVectorFunc::Function, problemDim::Int64)::Vector
    f::Vector = zeros(mesh.noOfNodes*problemDim)
    for elementNo ∈ 1:length(mesh.Elements[attribute])
        element::AbstractElement = mesh.Elements[attribute][elementNo]
        coordArray::Array{Float64,2} = RapidFEM.getCoordArray(mesh, element)
        shapeFunction::Array{ShapeFunction,1} = RapidFEM.feSpace!(FeSpace, element, mesh, RapidFEM.lagrange)
        f_local::Array{Float64,1} = localVectorFunc(problemfunction, problemDim, element, shapeFunction, coordArray)
        vNodes::Array{Int64} = RapidFEM.getVectorNodes(element, problemDim)
        f[vNodes] = f_local
    end
    return f
end

function getPermutionMatrix(vNodes::Array{Int64}, mesh::Mesh, problemDim::Int64)
    noOfVectorNodes::Int64 = mesh.noOfNodes*problemDim
    Pzeros::SparseMatrixCSC = spzeros(noOfVectorNodes, noOfVectorNodes)
    Pzeros[vNodes, vNodes] = sparse(I, length(vNodes), length(vNodes))
    P::SparseMatrixCSC = sparse(I, noOfVectorNodes, noOfVectorNodes)
    P = P - Pzeros
    println(Pzeros)
    return P, Pzeros
end

function applyDirichletBC!(A::SparseMatrixCSC, b::Vector, DirichletFunction::Function, attribute::Tuple{Int64, Int64}, mesh::Mesh, problemDim::Int64)
    for i ∈ 1:length(mesh.Elements[attribute...])
        element::AbstractElement = mesh.Elements[attribute...][i]
        vNodes::Array{Int64} = RapidFEM.getVectorNodes(element, problemDim)
        coordArray::Array{Float64} = RapidFEM.getCoordArray(mesh, element)
        for nodeNo ∈ 1:size(coordArray,2)
            b[vNodes[nodeNo:nodeNo+problemDim-1]] = DirichletFunction(coordArray[:,nodeNo])
        end
        P::SparseMatrixCSC, Pzeros::SparseMatrixCSC =  getPermutionMatrix(vNodes, mesh, problemDim)
        A[:] = P'*A*P
        A[:] = A + Pzeros
    end
end

function poissonEquation()
    mesh::Mesh = RapidFEM.readMesh("../test/Bar.msh")
    FeSpace = RapidFEM.createFeSpace()
    problemDim::Int64 = 1
    K::SparseMatrixCSC = assembleMatrix((3,4), FeSpace, mesh, local_lagrange_K, problemDim)
    source(x) = [0.0]
    f::Vector = assembleVector(source, (3,4), FeSpace, mesh, localSource, problemDim)
    neumann(x) = [1.0]
    f += assembleVector(neumann, (2,1), FeSpace, mesh, localNeumann, problemDim)
    DirichletFunction(x) = ones(1)
    applyDirichletBC!(K, f, DirichletFunction, (2,2), mesh, problemDim)
    return K,f
end
