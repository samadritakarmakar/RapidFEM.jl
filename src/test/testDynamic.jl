using LinearAlgebra, Plots
using RapidFEM, SparseArrays

include("testTimeInt.jl")
#=
M = spzeros(1,1)
M[1,1] = 1.0
C = spzeros(1,1)
C[1,1] = 2.0
K = spzeros(1,1)
K[1,1] = 3.0
u_n = zeros(1)
u_n[1] = 1.0
u̇_n = zeros(1)
u̇_n[1] = 1.0
SolutionArray = [u_n, u̇_n]
f_Array = [[1.0], [2.0], [4.0]]
Δt = 0.01
θ_Array = [0.5, 0.5]
for i ∈ 1:5
A_global, f_mean_global = SSpj_getFinal_A_b([M, C, K], SolutionArray, f_Array, Δt, θ_Array)
α_global = A_global\f_mean_global
#println(" A_global = " , A_global[1], " f_mean_global = ", f_mean_global[1], " α = ", α_global[1])
updateSolution!(SolutionArray, Δt, α_global)
println(SolutionArray, " ")
end=#

function testSspl()

    m, c, k = 1.0, 2.0, 3.0

    D, E, G = 6.0, 5.0, 4.0

    M = spzeros(1,1)
    M[1,1] = m
    C = spzeros(1,1)
    C[1,1] = c
    K = spzeros(1,1)
    K[1,1] = k

    usol(t) = D*t^2 + E*t + G
    u̇sol(t) = 2*D*t + E
    üsol(t) = 2*D

    u_n = usol(0.0)
    u̇_n = u̇sol(0.0)
    ü_n = üsol(0.0)
    SolutionArray = [[u_n], [u̇_n], [ü_n]]

    Δt = 0.05

    
    f_1 = -analyticalForceFunction(-0.05, m, c, k, D, E, G)
    f_2 = -analyticalForceFunction(0.0, m, c, k, D, E, G)
    #f_3 = analyticalForceFunction(0.05, m, c, k, D, E, G)
    f_Array = [[0], [f_1], [f_2]]
    println(f_Array)
    #first value of f_Array is anyway popped out by update_f! function

    θ_Array = [0.5, 0.5]

    for i ∈ 1:5
        f_current = [-analyticalForceFunction(i*Δt, m, c, k, D, E, G)]
        popfirst!(f_Array)
        push!(f_Array, f_current)
        println(f_Array)
        A_global, f_mean_global = SSpj_getFinal_A_b([M, C, K], SolutionArray, f_Array, Δt, θ_Array)
        α_global = A_global\f_mean_global
        updateSolution!(SolutionArray, Δt, α_global)
        println(SolutionArray, " ")
        println("t = ", i*Δt)
        println("Actual Solution: ", usol(i*Δt), " ", u̇sol(i*Δt), " ", üsol(i*Δt), "\n")
    end   
end

function get_hat_u_vals(u_n, u̇_n, ü_n, Δt::Real, β1::Real, β2::Real)
    ü_hat_n1 = -2.0/(β2*Δt^2)*u_n - 2.0/(β2*Δt)*u̇_n - (1.0 - β2)/(β2)*ü_n
    u̇_hat_n1 = -2.0*β1/(β2*Δt)*u_n + (1.0 - 2.0*β1/β2)*u̇_n + Δt*(1.0 - β1/β2)*Δt*ü_n
    return ü_hat_n1, u̇_hat_n1
end

#=function ImplicitNewmark(M_C_K::Vector, f, u_n, u̇_n, ü_n, Δt::Real, β1::Real, β2::Real)
    M = M_C_K[1]
    C = M_C_K[2]
    K = M_C_K[3]
    @assert β2 > 1e-8 "For Implicit Newmark, β2 must be greater than 1e-8"
    ü_hat_n1, u̇_hat_n1 = get_hat_u_vals(u_n, u̇_n, ü_n, Δt, β1, β2)
    A = 2.0/(β2*Δt^2)*M + 2.0*β1/(β2*Δt)*C + K
    b = f + M*(ü_hat_n1) + C*(u̇_hat_n1)
    u_n1 = -A\b
    ü_n1 = ü_hat_n1 + 2.0/(β2*Δt^2)*u_n1
    u̇_n1 = u̇_hat_n1 + 2.0*β1/(β2*Δt)*u_n1
    return u_n1, u̇_n1, ü_n1
end


function Newmark(M_C_K::Vector, f, u_n, u̇_n, ü_n, Δt::Real, β1::Real, β2::Real)
    M = M_C_K[1]
    C = M_C_K[2]
    K = M_C_K[3]
    
    u_breve_n1 = u_n + Δt*u̇_n + 0.5*(1.0 - β2)*Δt^2*ü_n
    u̇_breve_n1 = u̇_n + Δt*(1.0 - β1)*ü_n
    A = M + β1*Δt*C + 0.5*β2*Δt^2*K
    b = f + C*(u̇_breve_n1) + K*(u_breve_n1)
    ü_n1 = -A\b
    u̇_n1 = u̇_breve_n1 + β1*Δt*ü_n1
    u_n1 = u_breve_n1 + 0.5*Δt^2*β2*ü_n1
    return u_n1, u̇_n1, ü_n1
end=#


function testNewmark()
    m, c, k = 1.0, 2.0, 3.0
    D, E, G = 3.0, 2.0, 1.0

    M = m

    C = c

    K = k
    usol(t) = D*t^2 + E*t + G
    u̇sol(t) = 2*D*t + E
    üsol(t) = 2*D
    u_n = usol(0.0)
    u̇_n = u̇sol(0.0)
    ü_n = üsol(0.0)
    println("u_n = ", u_n, " u̇_n = ", u̇_n, " ü_n = ", ü_n)
    Δt = 0.05
    β1, β2 = (0.5, 0.5)
    for i ∈ 1:20
        f = -analyticalForceFunction(i*Δt, m, c, k, D, E, G)
        println("f = ", f)
        u_n1, u̇_n1, ü_n1 = ImplicitNewmark([M, C, K], f, u_n, u̇_n, ü_n,  Δt, β1, β2)
        #u_n1, u̇_n1, ü_n1 = Newmark([M, C, K], f, u_n, u̇_n, ü_n, Δt, β1, β2)
        println("t = ", i*Δt)
        println("u_n = ", u_n1, " u̇_n = ", u̇_n1, " ü_n = ", ü_n1)
        println("Actual Solution: ", usol(i*Δt), " ", u̇sol(i*Δt), " ", üsol(i*Δt), "\n")
        u_n, u̇_n, ü_n = u_n1, u̇_n1, ü_n1
    end
end


function NonLinearImplicitNewmark(u_n1, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f::T, u_n, u̇_n, ü_n, Δt::Real, β1::Real, β2::Real) where T
    @assert β2 > 1e-8 "For Implicit Newmark, β2 must be greater than 1e-8"
    ü_hat_n1, u̇_hat_n1 = get_hat_u_vals(u_n, u̇_n, ü_n, Δt, β1, β2)

    ü_n1 = ü_hat_n1 + 2.0/(β2*Δt^2)*u_n1
    u̇_n1 = u̇_hat_n1 + 2.0*β1/(β2*Δt)*u_n1

    if isa(f, Function)
        f = f(u_n1)
    end

    fTotal = f + M_vecfunc(ü_n1) + C_vecfunc(u̇_n1) + K_vecfunc(u_n1)

    return fTotal, u̇_n1,  ü_n1
end

function testImplicitNewmarkNonLinear()
    m, c, k = 1.0, 2.0, 3.0
    D, E, G = 3.0, 2.0, 1.0

    M = m

    C = c

    K = k
    usol(t) = D*t^2 + E*t + G
    u̇sol(t) = 2*D*t + E
    üsol(t) = 2*D
    u_n = usol(0.0)
    u̇_n = u̇sol(0.0)
    ü_n = üsol(0.0)

    M_vecfunc(ü_n) = M*ü_n
    C_vecfunc(u̇_n) = C*u̇_n
    K_vecfunc(u_n) = K*u_n

    Δt = 0.05

    β1, β2 = (0.5, 0.5)
    
    du = Inf
    u̇_n1, ü_n1 = 0.0, 0.0
    Δu = 0.0
    for i ∈ 1:10
        f = -analyticalForceFunction(i*Δt, m, c, k, D, E, G)
        f_func(u_n) = f
        println("f = ", f)
        fTotal = Inf
        u_n1 = u_n + Δu
        #finite difference
        iter = 0
        println("starting u_n1 = ", u_n1)
        while norm(fTotal) > 1e-8 || norm(du) > 1e-12
            fTotal, u̇_n1, ü_n1 = NonLinearImplicitNewmark(u_n1, M_vecfunc, C_vecfunc, K_vecfunc, f_func, u_n, u̇_n, ü_n, Δt, β1, β2)
            println("\nfTotal = ", fTotal)
            fTotalm1, u̇_n1_, ü_n1_ = NonLinearImplicitNewmark(u_n1 - 1e-4, M_vecfunc, C_vecfunc, K_vecfunc, f_func, u_n, u̇_n, ü_n,  Δt, β1, β2)
            
            fTotal1, u̇_n1_, ü_n1_ = NonLinearImplicitNewmark(u_n1 + 1e-4, M_vecfunc, C_vecfunc, K_vecfunc, f, u_n, u̇_n, ü_n, Δt, β1, β2)
            df = (fTotal1 - fTotalm1)/(2e-4)
            du = -df\fTotal
            println("du = ", du)
            u_n1 += du
            iter += 1
        end

        println("iter = ", iter)
        println("t = ", i*Δt)
        Δu = u_n1 - u_n
        u_n, u̇_n, ü_n = u_n1, u̇_n1, ü_n1
        println("u_n = ", u_n1, " u̇_n = ", u̇_n1, " ü_n = ", ü_n1)
        println("Actual Solution: ", usol(i*Δt), " ", u̇sol(i*Δt), " ", üsol(i*Δt), "\n")
        
    end

end

#=function getNewmark_Disp_Velocity(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    Δt::Real, β1::Real, β2::Real)

    u_n1 = u_n + Δt*u̇_n + 0.5*Δt^2*(1.0 - β2)*ü_n + 0.5*Δt^2*β2*ü_n1
    u̇_n1 = u̇_n + Δt*(1.0 - β1)*ü_n + β1*Δt*ü_n1
    return u_n1, u̇_n1
end

function NonLinearNewmark(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f::T, u_n::Union{Real, Vector}, 
    u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, Δt::Real, β1::Real, β2::Real) where T

    u_n1, u̇_n1 = getNewmark_Disp_Velocity(ü_n1, u_n, u̇_n, ü_n, Δt, β1, β2)
    
    #println("within Newmark: u_n1 = ", u_n1, " u̇_n1 = ", u̇_n1)

    if isa(f, Function)
        f = f(u_n1)
    end

    fTotal = f + M_vecfunc(ü_n1) + C_vecfunc(u̇_n1) + K_vecfunc(u_n1)

    return fTotal, u_n1, u̇_n1
end=#


function testNewmarkNonLinear()
    m, c, k = 1.0, 2.0, 3.0
    D, E, G = 3.0, 2.0, 1.0

    M = m

    C = c

    K = k
    usol(t) = D*t^2 + E*t + G
    u̇sol(t) = 2*D*t + E
    üsol(t) = 2*D
    u_n = usol(0.0)
    u̇_n = u̇sol(0.0)
    ü_n = üsol(0.0)

    M_vecfunc(ü_n) = M*ü_n
    C_vecfunc(u̇_n) = C*u̇_n
    K_vecfunc(u_n) = K*u_n

    Δt = 0.05

    β1 = 0.5
    β2 = 0.0
    dü = Inf
    Δü = 0.0
    u_n1, u̇_n1 = 0.0, 0.0
    u_nVal = Float64[]
    u̇_nVal = Float64[]
    ü_nVal = Float64[]
    t_nVal = Float64[]
    for i ∈ 1:5
        t = i*Δt
        f = -analyticalForceFunction(i*Δt, m, c, k, D, E, G)
        f_func(u_n, u̇_n, time) = f
        #println("u_n = $u_n u̇_n = $u̇_n ü_n = $ü_n")
        #println("f = ", f)
        fTotal = Inf
        ü_n1 = 0.0#ü_n + Δü
        #finite difference
        iter = 0
        #println("starting ü_n1 = ", ü_n1)
        while norm(fTotal) > 1e-8 || norm(dü) > 1e-12
            fTotal, u_n1, u̇_n1 = ImplicitNewmark(ü_n1, M_vecfunc, C_vecfunc, K_vecfunc, f_func, u_n, u̇_n, ü_n, t, t-Δt, β1, β2)
            #println("\nfTotal = ", fTotal)
            fTotalm1, u_n1_, u̇_n1_ = ImplicitNewmark(ü_n1 - 1e-4, M_vecfunc, C_vecfunc, K_vecfunc, f_func, u_n, u̇_n, ü_n,  t, t-Δt, β1, β2)
            fTotal1, u_n1_, u̇_n1_ = ImplicitNewmark(ü_n1 + 1e-4, M_vecfunc, C_vecfunc, K_vecfunc, f_func, u_n, u̇_n, ü_n, t, t-Δt, β1, β2)
            df = (fTotal1 - fTotalm1)/(2e-4)
            dü = -df\fTotal
            #println("dü = ", dü)
            ü_n1 += dü
            #println("during Newton: u_n = ", u_n1, " u̇_n = ", u̇_n1, " ü_n = ", ü_n1)
            
            iter += 1
        end
        push!(u_nVal, u_n1)
        push!(u̇_nVal, u̇_n1)
        push!(ü_nVal, ü_n1)
        push!(t_nVal, t)
        println("iter = ", iter)
        println("t = ", i*Δt)
        Δü = ü_n1 - ü_n
        u_n, u̇_n, ü_n = u_n1, u̇_n1, ü_n1
        println("u_n = ", u_n1, " u̇_n = ", u̇_n1, " ü_n = ", ü_n1)
        println("Actual Solution: ", usol(i*Δt), " ", u̇sol(i*Δt), " ", üsol(i*Δt), "\n")
    end
    println("ü_nVal = ", ü_nVal)
    accPlot = plot(label = "Acc", xlabel = "Time", ylabel = "Acceleration")
    plot!(accPlot, t_nVal, ü_nVal, label = "β1 = $β1, β2 = $β2")
end