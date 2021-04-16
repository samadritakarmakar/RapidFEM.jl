#=baseCoords =[0 0 0
0 1 0
0 0 1
1 0 0
 0 1/3 0
 0 2/3 0
 0 2/3 1/3
 0 1/3 2/3
 0 0 2/3
0 0 1/3
2/3 0 0
1/3 0 0
2/3 0 1/3
1/3 0 2/3
 2/3 1/3 0
 1/3 2/3 0
 0 1/3 1/3
 1/3 1/3 0
 1/3 0 1/3
 1/3 1/3 1/3]=#

# See wiki page:
# https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Conversion_between_barycentric_and_Cartesian_coordinates
function cartesian2barycentricTet(coords)
    cornerCoords::Array{Float64} = [0 0 0;
    0 1 0;
    0 0 1;
    1 0 0;]
    #x = [cornerCoords[1][1], cornerCoords[2][1], cornerCoords[3][1]]
    #y = [cornerCoords[1][2], cornerCoords[2][2], cornerCoords[3][2]]
    x = cornerCoords[:,1]
    y = cornerCoords[:,2]
    z = cornerCoords[:,3]
    #T = Matrix{Float64}(undef, 3, 3)
    T = zeros(4,4)
    T[1,1:4] = x[1:4]
    T[2,1:4] = y[1:4]
    T[3,1:4] = z[1:4]
    T[4,1:4] = ones(1, 4)
    #T[3,1:3] = ones(1, 3)
    baryCoords = Vector{Float64}(undef,4)
    baryCoords= T\[coords ;1]
    return baryCoords
end

function getTetBaryCentricCoordSet(baseCoords)
    baseCoordsBary::Array{Float64} = Array{Float64}(undef,size(baseCoords,1),4)
    i::UInt64 = 1
    for j ∈ 1: size(baseCoords,1)
        coord = baseCoords[j,:]
        coordPass = Vector{Float64}(undef, 3)
        coordPass[1] = coord[1]
        coordPass[2] = coord[2]
        coordPass[3] = coord[3]
        baseCoordsBary[j,:] = cartesian2barycentricTet(coordPass)
        #println(typeof(coordPass))
        #push!(baseCoordsBary[j],cartesian2barycentric(coordPass))
        i+=1
    end
    return (baseCoordsBary)
end

#getBaryCentricCoordSet(baseCoords)
