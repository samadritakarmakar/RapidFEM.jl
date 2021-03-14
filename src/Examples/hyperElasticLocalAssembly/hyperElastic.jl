function findStrain!(ϵ::Array{Float64, 1}, ∂ϕ_∂x::Array{Float64,2},
    solAtNodes::Array{Float64, 1}, problemDim::Int64)

    fill!(ϵ, 0.0)
    for a ∈ 1:size(∂ϕ_∂x, 1)
        for j ∈ 1:problemDim
            for i ∈ 1:problemDim
                ij::Int64 = LargeDeformations.getMandelIndex(i,j)
                ϵ[ij] += 0.5*(∂ϕ_∂x[a,j]*solAtNodes[problemDim*(a-1)+i]+∂ϕ_∂x[a,i]*solAtNodes[problemDim*(a-1)+j])
            end
        end
    end
end

function getCurrentCoordArray(coordArray::Array{Float64,2}, solAtNodes::Array{Float64, 1}, problemDim::Int64)
    currentCoordArray = zeros(size(coordArray)...)
    for a ∈ 1:size(coordArray, 2)
        for i ∈ 1:problemDim
            currentCoordArray[i,a] = coordArray[i, a] + solAtNodes[problemDim*(a-1)+i]
        end
    end
    return currentCoordArray
end

function get_∂u_∂X!(∂u_∂X::Array{Float64, 1}, solAtNodes::Array{Float64, 1}, ∂ϕ_∂X::Array{Float64,2}, problemDim::Int64)
    fill!(∂u_∂X, 0.0)
    for a ∈ 1:size(∂ϕ_∂X, 1)
        for J ∈ 1:problemDim
            @fastmath @simd for i ∈ 1:problemDim
                iJ::Int64 = LargeDeformations.getMandelIndex(i,J)
                ∂u_∂X[iJ] += ∂ϕ_∂X[a,J]*solAtNodes[problemDim*(a-1)+i]
            end
        end
    end
end

function get_∂u_∂X!(∂u_∂X::Array{Float64, 2}, solAtNodes::Array{Float64, 1}, ∂ϕ_∂X::Array{Float64,2}, problemDim::Int64)
    fill!(∂u_∂X, 0.0)
    for a ∈ 1:size(∂ϕ_∂X, 1)
        for J ∈ 1:problemDim
            @fastmath @simd for i ∈ 1:problemDim
                #iJ::Int64 = LargeDeformations.getMandelIndex(i,J)
                ∂u_∂X[i,J] += ∂ϕ_∂X[a,J]*solAtNodes[problemDim*(a-1)+i]
            end
        end
    end
end

function local_δE_S_Vector!(f::Vector, hyperElasticData::T, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T


    model::hyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]
    lastSoln::Array{Float64, 1} = hyperElasticData[3]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #ϵ::Array{Float64, 1} = zeros(model.ϵVoigtSize)
    ∂u_∂X::Array{Float64, 1} = zeros(problemDim^2)
    F = zeros(problemDim^2)
    Jacobian = 1.0
    E = zeros(problemDim^2)
    S = zeros(problemDim^2)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dX::Array{Float64,2} = ∂ξ_∂xFunc(∂X_∂ξ)
        dΩ::Float64 = dΩFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dX
        #findStrain!(ϵ, ∂ϕ_∂x,  solAtNodes, problemDim)
        get_∂u_∂X!(∂u_∂X, solAtNodes, ∂ϕ_∂X, problemDim)

        LargeDeformations.getDeformationGradient!(F, ∂u_∂X)
        #if ipNo == 1
        #    println("∂u_∂X = ", LargeDeformations.convert2DMandelToTensor(∂u_∂X))
        #end
        #println("E = ", LargeDeformations.getGreenStrain(F))
        #Jacobian = LargeDeformations.getJacobianDeformationGradient(F)
        LargeDeformations.getGreenStrain!(E, F)
        model.secondPiolaStress!(S, E, modelParams)
        #if ipNo == 1
        #    println("S = ", LargeDeformations.convert2DMandelToTensor(S))
        #end
        for a ∈ 1:noOfNodes
            for J::Int64 ∈ 1:problemDim
                ∂ϕ_∂X_a_J = ∂ϕ_∂X[a,J]
                for I::Int64 ∈ 1:problemDim
                    ∂ϕ_∂X_a_I = ∂ϕ_∂X[a,I]
                    IJ::Int64 = LargeDeformations.getMandelIndex(I, J)
                    c1::Float64 = 0.5#(i==j) ? 0.5 : 1.0
                    S_IJ =S[IJ]
                    @fastmath @simd for i ∈ 1:problemDim
                        iI::Int64 = LargeDeformations.getMandelIndex(i, I)
                        iJ::Int64 = LargeDeformations.getMandelIndex(i, J)
                        f[problemDim*(a-1)+i] += c1*(∂ϕ_∂X_a_I*F[iJ]+∂ϕ_∂X_a_J*F[iI])*
                        S_IJ*dΩ
                        #if  (i == 1 && S[IJ] > 39.0 && a==2 && f[problemDim*(a-1)+i] != 0.0)
                        #    println("∂ϕ_∂X[a,I]= ",∂ϕ_∂X[a,I], " F[iJ] = ", F[iJ], " ∂ϕ_∂X[a,J] = ", ∂ϕ_∂X[a,J], " F[iI] = ", F[iI], " S[IJ] = ", S[IJ])
                        #end
                    end
                end
            end
        end
    end
    return nothing
end

function local_∇v_σ_Vector!(f::Vector, hyperElasticData::T, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T


    model::hyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]
    lastSoln::Array{Float64, 1} = hyperElasticData[3]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, solAtNodes, problemDim)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    ϵ::Array{Float64, 1} = zeros(problemDim^2)
    ∂u_∂X::Array{Float64, 1} = zeros(problemDim^2)
    for ipNo ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(currentCoordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dX::Array{Float64,2} = ∂ξ_∂xFunc(∂X_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ::Float64 = dΩFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dX
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        #findStrain!(ϵ, ∂ϕ_∂x,  solAtNodes, problemDim)
        get_∂u_∂X!(∂u_∂X, solAtNodes, ∂ϕ_∂X, problemDim)
        F = LargeDeformations.getDeformationGradient(∂u_∂X)
        Jacobian = LargeDeformations.getJacobianDeformationGradient(F)
        σ = model.cauchyStress(F, modelParams)
        for a ∈ 1:noOfNodes
            for j::Int64 ∈ 1:problemDim
                for i::Int64 ∈ 1:problemDim
                    ij::Int64 = LargeDeformations.getMandelIndex(i, j)
                    c1::Float64 = 0.5#(i==j) ? 0.5 : 1.0
                    f[problemDim*(a-1)+i] += c1*∂ϕ_∂x[a,j]*σ[ij]*Jacobian*dΩ
                    f[problemDim*(a-1)+j] += c1*∂ϕ_∂x[a,i]*σ[ij]*Jacobian*dΩ
                end
            end
        end
    end
    #SmallStrainPlastic.updateStateDict!(stateDictBuffer, stateDictBufferCopy)
    return nothing
end

function local_δE_Cᵀ_ΔE!(𝕂::Array{Float64,2}, hyperElasticData::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T


    model::hyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]
    lastSoln::Array{Float64, 1} = hyperElasticData[3]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #ϵ::Array{Float64, 1} = zeros(model.ϵVoigtSize)
    ∂u_∂X::Array{Float64, 1} = zeros(problemDim^2)
    F = zeros(problemDim^2)
    Jacobian = 1.0
    E = zeros(problemDim^2)
    S = zeros(problemDim^2)
    ℂ = zeros(problemDim^2, problemDim^2)
    model.materialTangentTensor!(ℂ, rand(9), modelParams)
    ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[1].∂ϕ_∂ξ)
    ∂ξ_dX::Array{Float64,2} = ∂ξ_∂xFunc(∂X_∂ξ)
    dΩ::Float64 = dΩFunc(∂X_∂ξ, shapeFunction[1].ipData)
    ϕ::Array{Float64,1} = shapeFunction[1].ϕ
    #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
    ∂ϕ_∂X::Array{Float64} = shapeFunction[1].∂ϕ_∂ξ*∂ξ_dX
    for ipNo::Int64 ∈ 1:noOfIpPoints
        @inbounds ∂X_∂ξ = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        @inbounds ∂ξ_dX = ∂ξ_∂xFunc(∂X_∂ξ)
        dΩ = dΩFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        @inbounds ϕ = shapeFunction[ipNo].ϕ
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        @inbounds ∂ϕ_∂X = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dX
        #findStrain!(ϵ, ∂ϕ_∂x,  solAtNodes, problemDim)
        get_∂u_∂X!(∂u_∂X, solAtNodes, ∂ϕ_∂X, problemDim)
        LargeDeformations.getDeformationGradient!(F, ∂u_∂X)
        #println(F)
        #Jacobian = LargeDeformations.getJacobianDeformationGradient(F)
        LargeDeformations.getGreenStrain!(E, F)
        model.secondPiolaStress!(S, E, modelParams)
        #F_func(∂u_∂X) =
        #S_func(E_parm) = model.secondPiolaStress(E_parm, modelParams)
        #model.materialTangentTensor!(ℂ, E, modelParams)
        for b::Int64 ∈ 1:noOfNodes
            for a::Int64 ∈ 1:noOfNodes
                for L::Int64 ∈ 1:problemDim
                    ∂ϕ_∂X_b_L = ∂ϕ_∂X[b,L]
                    for K::Int64 ∈ 1:problemDim
                        ∂ϕ_∂X_b_K = ∂ϕ_∂X[b,K]
                        KL::Int64 = LargeDeformations.getMandelIndex(K, L)
                        #c2::Float64 = 0.5#(k==l) ? 0.5 : 1.0
                        for J::Int64 ∈ 1:problemDim
                            ∂ϕ_∂X_a_J = ∂ϕ_∂X[a,J]
                            ∂ϕ_∂X_b_J = ∂ϕ_∂X[b,J]
                            for I::Int64 ∈ 1:problemDim
                                IJ::Int64 = LargeDeformations.getMandelIndex(I, J)
                                #c1::Float64 = 0.5#(i==j) ? 0.5 : 1.0
                                ∂ϕ_∂X_a_I = ∂ϕ_∂X[a,I]
                                ∂ϕ_∂X_b_I = ∂ϕ_∂X[b,I]
                                S_IJ = S[IJ]
                                ℂ_IJKL = ℂ[IJ, KL]
                                @fastmath for j::Int64 ∈ 1:problemDim
                                    jK::Int64 = LargeDeformations.getMandelIndex(j, K)
                                    jL::Int64 = LargeDeformations.getMandelIndex(j, L)

                                    #𝕂[problemDim*(a-1)+j,problemDim*(b-1)+j] += (∂ϕ_∂X[a,I]*∂ϕ_∂X[b,J]+ ∂ϕ_∂X[a,J]*∂ϕ_∂X[b,I])*S[IJ]*dΩ
                                    𝕂[problemDim*(a-1)+j,problemDim*(b-1)+j] += (∂ϕ_∂X_a_I*∂ϕ_∂X_b_J + ∂ϕ_∂X_a_J*∂ϕ_∂X_b_I)*S_IJ*dΩ
                                    F_jL = F[jL]
                                    F_jK = F[jK]
                                    @fastmath for i::Int64 ∈ 1:problemDim
                                        iI::Int64 = LargeDeformations.getMandelIndex(i, I)
                                        iJ::Int64 = LargeDeformations.getMandelIndex(i, J)

                                        #=𝕂[problemDim*(a-1)+i,problemDim*(b-1)+j] += 0.25*
                                        (∂ϕ_∂X[a,I]*F[iJ]+∂ϕ_∂X[a,J]*F[iI])*
                                        ℂ[IJ, KL]*
                                        (∂ϕ_∂X[b,K]*F[jL]+∂ϕ_∂X[b,L]*F[jK])*dΩ=#
                                        𝕂[problemDim*(a-1)+i,problemDim*(b-1)+j] += 0.25*
                                        (∂ϕ_∂X_a_I*F[iJ]+∂ϕ_∂X_a_J*F[iI])*
                                        ℂ_IJKL*
                                        (∂ϕ_∂X_b_K*F[jL]+∂ϕ_∂X_b_L*F[jK])*dΩ
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    #SmallStrainPlastic.updateStateDict!(stateDictBuffer, stateDictBufferCopy)
    return nothing
end

function local_∇v_Cᵀ_∇u!(K::Array{Float64,2}, hyperElasticData::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...)where T

    model::hyperElasticModel = hyperElasticData[1]
    modelParams::Tuple  = hyperElasticData[2]
    lastSoln::Array{Float64, 1} = hyperElasticData[3]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, solAtNodes, problemDim)
    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #ϵ::Array{Float64, 1} = zeros(model.ϵVoigtSize)
    ∂u_∂X::Array{Float64, 1} = zeros(problemDim, problemDim)
    for ipNo::Int64 ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(currentCoordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dX::Array{Float64,2} = ∂ξ_∂xFunc(∂X_∂ξ)
        ∂ξ_dx::Array{Float64,2} = ∂ξ_∂xFunc(∂x_∂ξ)
        dΩ::Float64 = dΩFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        #X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        ∂ϕ_∂X::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dX
        ∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        #findStrain!(ϵ, ∂ϕ_∂x,  solAtNodes, problemDim)
        get_∂u_∂X!(∂u_∂X, solAtNodes, ∂ϕ_∂X, problemDim)
        F = LargeDeformations.getDeformationGradient(∂u_∂X)
        Jacobian = LargeDeformations.getJacobianDeformationGradient(F)
        𝕔 = model.spatialTangentTensor(F, modelParams)
        for b::Int64 ∈ 1:noOfNodes
            for a::Int64 ∈ 1:noOfNodes
                for l::Int64 ∈ 1:problemDim
                    for k::Int64 ∈ 1:problemDim
                        kl::Int64 = LargeDeformations.getMandelIndex(k, l)
                        c2::Float64 = 0.5#(k==l) ? 0.5 : 1.0
                        #c2 = 0.5
                        for j::Int64 ∈ 1:problemDim
                            for i::Int64 ∈ 1:problemDim
                                ij::Int64 = LargeDeformations.getMandelIndex(i, j)
                                c1::Float64 = 0.5#(i==j) ? 0.5 : 1.0
                                #c1 = 0.5
                                K[problemDim*(a-1)+i,problemDim*(b-1)+k] += c1*c2*∂ϕ_∂x[a,j]*𝕔[ij,kl]*∂ϕ_∂x[b,l]*Jacobian*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+l] += c1*c2*∂ϕ_∂x[a,i]*𝕔[ij,kl]*∂ϕ_∂x[b,k]*Jacobian*dΩ
                                K[problemDim*(a-1)+j,problemDim*(b-1)+k] += c1*c2*∂ϕ_∂x[a,i]*𝕔[ij,kl]*∂ϕ_∂x[b,l]*Jacobian*dΩ
                                K[problemDim*(a-1)+i,problemDim*(b-1)+l] += c1*c2*∂ϕ_∂x[a,j]*𝕔[ij,kl]*∂ϕ_∂x[b,k]*Jacobian*dΩ
                            end
                        end
                    end
                end
            end
        end
    end
    #SmallStrainPlastic.updateStateDict!(stateDictBuffer, stateDictBufferCopy)
    return nothing
end

function localReferenceSource!(S::Vector, hyperElasticData::T, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T


    sourceFunc = hyperElasticData[1]
    lastSoln::Array{Float64, 1} = hyperElasticData[2]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, solAtNodes, problemDim)

    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #S::Vector = zeros(noOfNodes*problemDim)
    ∂u_∂X::Array{Float64, 1} = zeros(problemDim^2)
    for ipNo ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dX::Array{Float64,2} = ∂ξ_∂xFunc(∂X_∂ξ)
        dΩ::Float64 = dΩFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        ∂ϕ_∂X::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dX
        get_∂u_∂X!(∂u_∂X, solAtNodes, ∂ϕ_∂X, problemDim)
        F = LargeDeformations.getDeformationGradient(∂u_∂X)
        J = LargeDeformations.getJacobianDeformationGradient(F)
        X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        s::Array{Float64,1} = sourceFunc(X; kwargs4function...)
        for a ∈ 1:noOfNodes
            @fastmath @simd for i ∈ 1:problemDim
                S[problemDim*(a-1)+i] += ϕ[a]*s[i]*dΩ
            end
        end
    end
    return nothing
end

function localCurrentSource!(S::Vector, hyperElasticData::T, problemDim::Int64,
    element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T


    sourceFunc = hyperElasticData[1]
    lastSoln::Array{Float64, 1} = hyperElasticData[2]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, solAtNodes, problemDim)

    ∂ξ_∂xFunc::Function = getFunction_∂ξ_∂x(element)
    dΩFunc::Function = getFunction_dΩ(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #S::Vector = zeros(noOfNodes*problemDim)
    ∂u_∂X::Array{Float64, 1} = zeros(problemDim^2)
    for ipNo ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂ξ_dX::Array{Float64,2} = ∂ξ_∂xFunc(∂X_∂ξ)
        dΩ::Float64 = dΩFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        ∂ϕ_∂X::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dX
        get_∂u_∂X!(∂u_∂X, solAtNodes, ∂ϕ_∂X, problemDim)
        F = LargeDeformations.getDeformationGradient(∂u_∂X)
        J = LargeDeformations.getJacobianDeformationGradient(F)
        x::Array{Float64,1} = getInterpolated_x(currentCoordArray, ϕ)
        s::Array{Float64,1} = sourceFunc(x; kwargs4function...)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                S[problemDim*(a-1)+i] += ϕ[a]*s[i]*J*dΩ
            end
        end
    end
    return nothing
end

function localReferenceNeumann!(Nm::Vector, hyperElasticData::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    neumannFunc = hyperElasticData[1]
    lastSoln::Array{Float64, 1} = hyperElasticData[2]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, solAtNodes, problemDim)

    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #Nm::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        dS::Float64 = dSFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        nm::Array{Float64,1} = neumannFunc(X; kwargs4function...)
        for a ∈ 1:noOfNodes
            @fastmath @simd for i ∈ 1:problemDim
                Nm[problemDim*(a-1)+i] += ϕ[a]*nm[i]*dS
            end
        end
    end
    return nothing
end

function localCurrentNeumann!(Nm::Vector, hyperElasticData::T,
    problemDim::Int64, element::AbstractElement, elementNo::Int64, shapeFunction::Array{ShapeFunction},
    coordArray::Array{Float64,2}; kwargs4function...) where T

    neumannFunc = hyperElasticData[1]
    lastSoln::Array{Float64, 1} = hyperElasticData[2]

    solAtNodes::Array{Float64, 1} = getSolAtElement(lastSoln, element, problemDim)
    currentCoordArray = getCurrentCoordArray(coordArray, solAtNodes, problemDim)

    dSFunc::Function = getFunction_dS(element)
    noOfIpPoints::Int64 = length(shapeFunction)
    noOfNodes::Int64 = size(shapeFunction[1].∂ϕ_∂ξ,1)
    #Nm::Vector = zeros(noOfNodes*problemDim)
    for ipNo ∈ 1:noOfIpPoints
        ∂X_∂ξ::Array{Float64,2} = get_∂x_∂ξ(coordArray, shapeFunction[ipNo].∂ϕ_∂ξ)
        ∂x_∂ξ::Array{Float64,2} = get_∂x_∂ξ(currentCoordArray, shapeFunction[ipNo].∂ϕ_∂ξ)

        ds::Float64 = dSFunc(∂x_∂ξ, shapeFunction[ipNo].ipData)
        dS::Float64 = dSFunc(∂X_∂ξ, shapeFunction[ipNo].ipData)
        ϕ::Array{Float64,1} = shapeFunction[ipNo].ϕ
        #∂ϕ_∂x::Array{Float64} = shapeFunction[ipNo].∂ϕ_∂ξ*∂ξ_dx
        #get_∂u_∂X!(∂u_∂X, solAtNodes, ∂ϕ_∂X, problemDim)
        #F = LargeDeformations.getDeformationGradient(∂u_∂X)
        #J = LargeDeformations.getJacobianDeformationGradient(F)
        x::Array{Float64,1} = getInterpolated_x(currentCoordArray, ϕ)
        X::Array{Float64,1} = getInterpolated_x(coordArray, ϕ)
        nm::Array{Float64,1} = neumannFunc(x; kwargs4function...)
        for a ∈ 1:noOfNodes
            for i ∈ 1:problemDim
                Nm[problemDim*(a-1)+i] += ϕ[a]*nm[i]*ds
            end
        end
    end
    return nothing
end
