"""For non-linear problems, this function returns the next set of possible displacement, velocity, given the current acceleration guess in 
a Newton-Raphson iteration."""
function getNewmark_Disp_Velocity(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    Δt::Real, β1::Real, β2::Real)

    u_n1 = u_n + Δt*u̇_n + 0.5*Δt^2*(1.0 - β2)*ü_n + 0.5*Δt^2*β2*ü_n1
    u̇_n1 = u̇_n + Δt*(1.0 - β1)*ü_n + β1*Δt*ü_n1
    return u_n1, u̇_n1
end

"""This function returns the total force acting on the system, given the current acceleration guess in a Newton-Raphson iteration.
The mass, damping, and stiffness matrices are represented by the functions M_vecfunc(ü_n1), C_vecfunc(u̇_n1), and K_vecfunc(u_n1), respectively. The
external force acting on the system is represented by f. The current displacement, velocity, and acceleration are represented by u_n,
u̇_n, and ü_n, respectively. The time step is represented by Δt. The Newmark parameters are represented by β1 and β2."""
function NonLinearNewmark(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f::T, u_n::Union{Real, Vector}, 
    u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, Δt::Real, β1::Real, β2::Real) where T

    u_n1, u̇_n1 = getNewmark_Disp_Velocity(ü_n1, u_n, u̇_n, ü_n, Δt, β1, β2)
    
    #println("within Newmark: u_n1 = ", u_n1, " u̇_n1 = ", u̇_n1)

    if isa(f, Function)
        f = f(u_n1)
    end

    fTotal = f + M_vecfunc(ü_n1) + C_vecfunc(u̇_n1) + K_vecfunc(u_n1)

    return fTotal, u_n1, u̇_n1
end