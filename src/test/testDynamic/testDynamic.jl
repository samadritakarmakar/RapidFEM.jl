using LinearAlgebra, Plots, LaTeXStrings
using RapidFEM, SparseArrays
include("testTimeInt.jl")

#Plots.pgfplotsx()

function testNewmarkNonLinear()
    m, c, k = 1.0, 2.0, 3.0


    tTotal = 10.0
    

    model::String = "oscillator"
    wavelength = nothing
    if model == "polyO2"
        D, E, G = 3.0, 2.0, 1.0
        paramsPoly = (m, c, k, D, E, G)
        usolFunc, u̇solFunc, üsolFunc, analyticalF_Function = polynomialModel2ndOrder
        params = paramsPoly
    elseif model == "oscillator"
        D = 1.0
        paramsOscillator = (m, c, k, D)
        usolFunc, u̇solFunc, üsolFunc, analyticalF_Function = oscillatorModel
        params = paramsOscillator
        wavelength = 2*π/D
    elseif model == "rampFlat"
        u_total, TaccEnd, Ttotal = 10.0, 0.25*tTotal, tTotal
        paramsRampFlat = (m, c, k, u_total, TaccEnd, Ttotal)
        usolFunc, u̇solFunc, üsolFunc, analyticalF_Function = rampFlatModel
        params = paramsRampFlat
    end



    usol(t) = usolFunc(t, params)
    u̇sol(t) = u̇solFunc(t, params)
    üsol(t) = üsolFunc(t, params)
    analyticalForceFunction(t) = analyticalF_Function(t, params)
    u_n = usol(0.0)
    u̇_n = u̇sol(0.0)
    ü_n = üsol(0.0)

    M_vecfunc(ü_n) = m*ü_n
    C_vecfunc(u̇_n) = c*u̇_n
    K_vecfunc(u_n) = k*u_n

    Δt = 1.2

    β1 = 0.5
    β2 = 0.5
    dü = Inf
    Δü = 0.0
    u_n1, u̇_n1 = 0.0, 0.0
    tVec = [0.0]
    uSolVec = [u_n]
    u̇SolVec = [u̇_n]
    üSolVec = [ü_n]

    uAnalytical = [usol(0.0)]
    u̇Analytical = [u̇sol(0.0)]
    üAnalytical = [üsol(0.0)]

    i = 1.0
    tEnd = 0.0
    while i*Δt <= tTotal
        f = -analyticalForceFunction(i*Δt)
        push!(tVec, i*Δt)
        f_func(u_n, u̇_n, t) = f
        #println("u_n = $u_n u̇_n = $u̇_n ü_n = $ü_n")
        #println("f = ", f)
        fTotal = Inf
        ü_n1 = 0.0#ü_n + Δü
        #finite difference
        iter = 0
        #println("starting ü_n1 = ", ü_n1)
        while norm(fTotal) > 1e-8 || norm(dü) > 1e-8
            fTotal, u_n1, u̇_n1 = ImplicitNewmark(ü_n1, M_vecfunc, C_vecfunc, K_vecfunc, f_func, u_n, u̇_n, ü_n, i*Δt, (i-1)*Δt, β1, β2)
            #println("\nfTotal = ", fTotal)
            fTotalm1, u_n1_, u̇_n1_ = ImplicitNewmark(ü_n1 - 1e-4, M_vecfunc, C_vecfunc, K_vecfunc, f_func, u_n, u̇_n, ü_n,  i*Δt, (i-1)*Δt, β1, β2)
            fTotal1, u_n1_, u̇_n1_ = ImplicitNewmark(ü_n1 + 1e-4, M_vecfunc, C_vecfunc, K_vecfunc, f, u_n, u̇_n, ü_n, i*Δt, (i-1)*Δt, β1, β2)
            df = (fTotal1 - fTotalm1)/(2e-4)
            dü = -df\fTotal
            println("dü = ", dü)
            ü_n1 += dü
            #println("during Newton: u_n = ", u_n1, " u̇_n = ", u̇_n1, " ü_n = ", ü_n1)
            
            iter += 1
        end
        println("iter = ", iter)
        println("t = ", i*Δt)
        Δü = ü_n1 - ü_n
        u_n, u̇_n, ü_n = u_n1, u̇_n1, ü_n1
        push!(uSolVec, u_n)
        push!(u̇SolVec, u̇_n)
        push!(üSolVec, ü_n)

        println("u_n = ", u_n1, " u̇_n = ", u̇_n1, " ü_n = ", ü_n1)
        println("Actual Solution: ", usol(i*Δt), " ", u̇sol(i*Δt), " ", üsol(i*Δt), "\n")
        push!(uAnalytical, usol(i*Δt))
        push!(u̇Analytical, u̇sol(i*Δt))
        push!(üAnalytical, üsol(i*Δt))
        tEnd = i*Δt
        i += 1.0
    end
    println("tEnd = ", tEnd)
    println("wavelength = ", wavelength)
    if model == "oscillator"
        println("Soln at wavelength= " , usol(2/D))
    end
    analyticTime = collect(LinRange(0.0, tEnd, 1000))

    #create directory for storing plots as per Δt
    dirName = "newmarkPlots-$(Δt)"
    if !isdir(dirName)
        mkdir(dirName)

    end

    plotDisp = plot(xlabel = "Time Δt=$Δt", ylabel = L"u", legend = :topleft, framestyle = :box,
    xtickfontsize=10,ytickfontsize=10, xguidefontsize=10, yguidefontsize=10,legendfontsize=10, legendfonthalign = :left,
    linewidth=1, thickness_scaling = 1)
    plot!(plotDisp, analyticTime, usol.(analyticTime), label = "Analytic")
    plot!(plotDisp, tVec, uSolVec, label = L"Numeric: $\beta_2=$%$β2", line = :dash, marker = :circle)
    savefig(plotDisp, "$dirName/newmarkDisp-Model-$(model)_Beta1-$(β1)_Beta2-$(β2)_DelT-$(Δt).svg")

    plotDispError = plot(xlabel = "Time Δt=$Δt", ylabel = L"u - u_{Anlyt}", legend = :topright, framestyle = :box,
    xtickfontsize=10,ytickfontsize=10, xguidefontsize=10, yguidefontsize=10,legendfontsize=10, legendfonthalign = :left,
    linewidth=1, thickness_scaling = 1)
    plot!(plotDispError, tVec, uSolVec .- uAnalytical, label = L"Error: $\beta_2=$%$β2")
    savefig(plotDispError, "$dirName/newmarkDispError-Model-$(model)_Beta1-$(β1)_Beta2-$(β2)_DelT-$(Δt).svg")

    plotVel = plot(xlabel = "Time Δt=$Δt", ylabel = L"\dot{u}", legend = :bottomright, framestyle = :box,
    xtickfontsize=10,ytickfontsize=10, xguidefontsize=10, yguidefontsize=10,legendfontsize=10, legendfonthalign = :left,
    linewidth=1, thickness_scaling = 1)
    plot!(plotVel, analyticTime, u̇sol.(analyticTime), label = "Analytic")
    plot!(plotVel, tVec, u̇SolVec, label = L"Numeric: $\beta_2=$%$β2", line = :dash, marker = :circle)
    savefig(plotVel, "$dirName/newmarkVel-Model-$(model)_Beta1-$(β1)_Beta2-$(β2)_DelT-$(Δt).svg")

    plotVelError = plot(xlabel = "Time Δt=$Δt", ylabel = L"\dot{u} - \dot{u}_{Anlyt}", legend = :bottomright, framestyle = :box,
    xtickfontsize=10,ytickfontsize=10, xguidefontsize=10, yguidefontsize=10,legendfontsize=10, legendfonthalign = :left,
    linewidth=1, thickness_scaling = 1)
    plot!(plotVelError, tVec, u̇SolVec .- u̇Analytical, label = L"Error: $\beta_2=$%$β2")
    savefig(plotVelError, "$dirName/newmarkVelError-Model-$(model)_Beta1-$(β1)_Beta2-$(β2)_DelT-$(Δt).svg")

    plotAcc = plot(xlabel = "Time Δt=$Δt", ylabel = L"\ddot{u}", legend = :bottomright, framestyle = :box,
    xtickfontsize=10,ytickfontsize=10, xguidefontsize=10, yguidefontsize=10,legendfontsize=10, legendfonthalign = :left,
    linewidth=1, thickness_scaling = 1)
    plot!(plotAcc, analyticTime, üsol.(analyticTime), label = "Analytic")
    plot!(plotAcc, tVec, üSolVec, label = L"Numeric Soln: $\beta_2=$%$β2", line = :dash, marker = :circle)
    savefig(plotAcc, "$dirName/newmarkAcc-Model-$(model)_Beta1-$(β1)_Beta2-$(β2)_DelT-$(Δt).svg")

    plotAccError = plot(xlabel = "Time Δt=$Δt", ylabel = L"\ddot{u} - \ddot{u}_{Anlyt}", legend = :bottomright, framestyle = :box,
    xtickfontsize=10,ytickfontsize=10, xguidefontsize=10, yguidefontsize=10,legendfontsize=10, legendfonthalign = :left,
    linewidth=1, thickness_scaling = 1)
    plot!(plotAccError, tVec, üSolVec .- üAnalytical, label = L"Error: $\beta_2=$%$β2")
    savefig(plotAccError, "$dirName/newmarkAccError-Model-$(model)_Beta1-$(β1)_Beta2-$(β2)_DelT-$(Δt).svg")
end
