function simpleNLsolve(assemble_f::Function, assemble_J::Function, initSoln::AbstractVector;
    xtol = 1e-8, ftol = 1e-8, iterations = 30)
    f = assemble_f(initSoln)
    ΔSoln = zeros(length(f))
    iter = 0
    while (norm(f)> ftol || norm(ΔSoln) > xtol) && iter < iterations
        J = assemble_J(initSoln)
        ΔSoln = J\f
        initSoln -= ΔSoln
        f = assemble_f(initSoln)
        iter += 1
        #println("norm(ΔSoln) = ", norm(ΔSoln))
    end
    if (iter>=iterations)
        @warn "simpleNLsolve exited without convergence!"
    end
    return initSoln
end

function simpleNLsolve(assemble_fJ::Function, initSoln::AbstractVector;
    xtol = 1e-8, ftol = 1e-8, iterations = 30)
    f, J = assemble_fJ(initSoln)
    ΔSoln = zeros(length(f))
    iter = 0
    while (norm(f)> ftol || norm(ΔSoln) > xtol) && iter < iterations
        ΔSoln = J\f
        initSoln -= ΔSoln
        f, J = assemble_fJ(initSoln)
        iter += 1
        #println("norm(ΔSoln) = ", norm(ΔSoln))
    end
    if (iter>=iterations)
        @warn "simpleNLsolve exited without convergence!"
    end
    return initSoln
end
