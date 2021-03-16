struct Convergence
    f_norm::Array{Float64, 1}
    Δx_norm::Array{Float64, 1}
    relNorm::Array{Float64, 1}
end


function simpleNLsolve(assemble_f::Function, assemble_J::Function, initSoln::AbstractVector;
    xtol = 1e-8, ftol = 1e-8, relTol= 1e-8, iterations = 30, skipJacobian = 1, printConvergence = false)
    f = assemble_f(initSoln)
    J = assemble_J(initSoln)
    ΔSoln = zeros(length(f))
    ΔSoln .= J\f
    initSoln -= ΔSoln
    iter1ResdualNorm = norm(initSoln)
    relNorm = 1.0
    f_normArray::Array{Float64, 1} = zeros(0)
    Δx_normArray::Array{Float64, 1} = zeros(0)
    relNormArray::Array{Float64, 1} = zeros(0)
    iter = 0
    while ((norm(f)> ftol || norm(ΔSoln) > xtol) && relNorm > relTol) && iter < iterations
        iter += 1
        f = assemble_f(initSoln)
        if (mod(iter,skipJacobian) == 0)
            J = assemble_J(initSoln)
        end
        ΔSoln .= J\f
        initSoln -= ΔSoln
        if printConvergence
            norm_f = norm(f)
            Δx_norm =norm(ΔSoln)
             relNorm = norm_f/iter1ResdualNorm
            println("\nnorm(f) = ", norm_f, " norm(Δx) = ", norm(ΔSoln), " Iteration = ", iter, " relNorm = ", norm_f/iter1ResdualNorm ,"\n")
            push!(f_normArray, norm_f)
            push!(Δx_normArray, Δx_norm)
            push!(relNormArray, relNorm)
        end
        #println("norm(ΔSoln) = ", norm(ΔSoln))
    end
    if (iter>=iterations)
        @warn "simpleNLsolve exited without convergence!"
    end

    return initSoln, Convergence(f_normArray, Δx_normArray, relNormArray)
end

function simpleNLsolve(assemble_fJ::Function, initSoln::AbstractVector;
    xtol = 1e-8, ftol = 1e-8, iterations = 30, printConvergence = false)
    f, J = assemble_fJ(initSoln)
    ΔSoln = zeros(length(f))
    iter = 0
    while (norm(f)> ftol || norm(ΔSoln) > xtol) && iter < iterations
        ΔSoln .= J\f
        initSoln -= ΔSoln
        f, J = assemble_fJ(initSoln)
        iter += 1
        if printConvergence
            println("\nnorm(f) = ", norm(f), " norm(Δx) = ", norm(ΔSoln), " Iteration = ", iter, "\n")
        end
        #println("norm(ΔSoln) = ", norm(ΔSoln))
    end
    if (iter>=iterations)
        @warn "simpleNLsolve exited without convergence!"
    end
    return initSoln
end
