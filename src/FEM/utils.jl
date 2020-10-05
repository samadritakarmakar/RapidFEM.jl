"""This function finds the size of the side of the element which has the lowest inclination
to the given velocity vector

    elmntSizeAlongVel(mesh, element, velocity)
"""
function elmntSizeAlongVel(mesh::Mesh, element::AbstractElement, velocity::Array{Float64,1})
    coordArrayTrans::Array{Float64,2} = getCoordArray(mesh, element)'
    noOfRows = size(coordArrayTrans, 1)
    eye::Array{Float64, 2} = Array{Float64, 2}(I, noOfRows, noOfRows)
    permMatrix::Array{Float64, 2} = zeros(noOfRows, noOfRows)
    permMatrix[1,:] = eye[3,:]
    for i ∈ 2:noOfRows
        permMatrix[i,:] = eye[i-1,:]
    end
    coordArrayPerm::Array{Float64,2} = permMatrix*coordArrayTrans
    vecMat::Array{Float64, 2} = coordArrayTrans-coordArrayPerm
    vecNorm::Array{Float64, 1} = Array{Float64, 1}(undef, size(vecMat, 1))
    for i ∈ size(vecMat,1)
        vecNorm[i] = norm(vecMat[i,:])
    end
    cosθ::Array{Float64, 1} = vecMat*velocity/(vecNorm*norm(velocity))
    maxCosθ::Float64, n::Int64 = findmax(cosθ)
    return coordArrayTrans[n, :]'
end
