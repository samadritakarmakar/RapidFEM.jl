function getSolutionVectorIn3d(solution::Vector{Float64}, activeDims::Vector{Int64})
    actualDim = sum(activeDims)
    noOfNodes = Int64(length(solution)/actualDim)
    newSolution = zeros(noOfNodes*3)
    Threads.@threads for nodeNo ∈ 1:noOfNodes
        currentDim3d = 1
        currentDimXd = 1
        for activeDim ∈ activeDims
            if activeDim != 0
                newSolution[3*(nodeNo-1)+currentDim3d] = solution[actualDim*(nodeNo -1)+currentDimXd]
                currentDimXd += 1
            end
            currentDim3d += 1
        end
    end
    return newSolution
end