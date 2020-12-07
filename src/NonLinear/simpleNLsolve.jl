function simpleNLsolve(assemble_f::Function, assemble_J::Function, initSoln::AbstractVector;
    xtol = 1e-8, ftol = 1e-8, iterations = 30, skipJacobian = 1, printConvergence = false)
    f = assemble_f(initSoln)
    J = assemble_J(initSoln)
    ΔSoln = zeros(length(f))
    ΔSoln = J\f
    initSoln -= ΔSoln
    iter = 0
    while (norm(f)> ftol || norm(ΔSoln) > xtol) && iter < iterations
        iter += 1
        f = assemble_f(initSoln)
        if (mod(iter,skipJacobian) == 0)
            J = assemble_J(initSoln)
        end
        ΔSoln = J\f
        initSoln -= ΔSoln
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

function simpleNLsolve(assemble_fJ::Function, initSoln::AbstractVector;
    xtol = 1e-8, ftol = 1e-8, iterations = 30, printConvergence = false)
    f, J = assemble_fJ(initSoln)
    ΔSoln = zeros(length(f))
    iter = 0
    while (norm(f)> ftol || norm(ΔSoln) > xtol) && iter < iterations
        ΔSoln = J\f
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
