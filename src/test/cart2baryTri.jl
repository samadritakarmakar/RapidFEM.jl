#=baseCoords =
[0.0 0.0 0.0;
1.0 0.0 0.0;
0.0 1.0 0.0;
1.0/3.0 0.0 0.0;
2.0/3.0 0.0 0.0;
2.0/3.0 1.0/3.0 0.0;
1.0/3.0 2.0/3.0 0.0;
0.0 2.0/3.0 0.0;
0.0 1.0/3.0 0.0;
1.0/3.0 1.0/3.0 0.0]=#

# See wiki page:
# https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Conversion_between_barycentric_and_Cartesian_coordinates
function cartesian2barycentricTri(coords)
    cornerCoords::Array{Float64} = [0 0 0; 1 0 0; 0 1 0]
    #x = [cornerCoords[1][1], cornerCoords[2][1], cornerCoords[3][1]]
    #y = [cornerCoords[1][2], cornerCoords[2][2], cornerCoords[3][2]]
    x = cornerCoords[:,1]
    y = cornerCoords[:,2]
    T = Matrix{Float64}(undef, 3, 3)
    T[1,1:3] = x[1:3]
    T[2,1:3] = y[1:3]
    T[3,1:3] = ones(1, 3)
    baryCoords = Vector{Float64}(undef,3)
    baryCoords= T\[coords ;1]
    return baryCoords
end

function getTriBaryCentricCoordSet(baseCoords::Array{Float64})
    baseCoordsBary::Array{Float64} = Array{Float64}(undef,size(baseCoords,1),3)
    i::UInt64 = 1
    for j ∈ 1: size(baseCoords,1)
        coord = baseCoords[j,:]
        coordPass = Vector{Float64}(undef, 2)
        coordPass[1] = coord[1]
        coordPass[2] = coord[2]
        baseCoordsBary[j,:] = cartesian2barycentricTri(coordPass)
        #println(typeof(coordPass))
        #push!(baseCoordsBary[j],cartesian2barycentric(coordPass))
        i+=1
    end
    return (baseCoordsBary)
end

#getBaryCentricCoordSet(baseCoords)
