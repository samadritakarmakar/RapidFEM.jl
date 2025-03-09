"""For non-linear problems, this function returns the next set of possible displacement, velocity, given the current acceleration guess in 
a Newton-Raphson iteration. WARNING: β2 parameter is here is 2*β parameter in the literature. β1 = γ as in the literature."""
function getNewmark_Disp_Velocity(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    Δt::Real, β1::Real, β2::Real)

    u_n1 = u_n + Δt*u̇_n + 0.5*Δt^2*(1.0 - β2)*ü_n + 0.5*Δt^2*β2*ü_n1
    u̇_n1 = u̇_n + Δt*(1.0 - β1)*ü_n + β1*Δt*ü_n1
    return u_n1, u̇_n1
end

"""Given the current displacement: u_n1, velocity: u̇_n1, and the previous displacement: u_n, velocity: u̇_n, this function returns the
accelerations at the previous: ü_n, and current time steps: ü_n1. The time step is represented by Δt. The Newmark parameters are represented
by β1 and β2. WARNING: β2 parameter is here is 2*β parameter in the literature. β1 = γ as in the literature."""
function getNewmarkAccsFromDispVel(u_n1::Union{Real, Vector}, u̇_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, Δt::Real, β1::Real, β2::Real)
    if β1 == β2
        error("β1 and β2 cannot be equal to find accelerations with this function.")
    end
    ü_n = ((2.0*β1*(u_n1 - u_n)/Δt - u̇_n) - β2*u̇_n1 + β2*u̇_n)/((β1 - β2)*Δt)
    ü_n1 = ((u̇_n1 - u̇_n)/Δt - ü_n)/β1 + ü_n
    return ü_n, ü_n1
end


"""This function returns the total force acting on the system, given the current acceleration guess in a Newton-Raphson iteration.
The mass, damping, and stiffness matrices are represented by the functions M_vecfunc(ü_n1), C_vecfunc(u̇_n1), and K_vecfunc(u_n1), respectively. The
external force acting on the system is represented by f. The current displacement, velocity, and acceleration are represented by u_n,
u̇_n, and ü_n, respectively. The time is represented by t_n and t_n1 (t_n+1). The Newmark parameters are represented by β1 and β2.
WARNING: β2 parameter is here is 2*β parameter in the literature. β1 = γ as in the literature."""
function ImplicitNewmark(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f::T, u_n::Union{Real, Vector}, 
    u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, t_n1::Real, t_n::Real, β1::Real, β2::Real) where T

    u_n1, u̇_n1 = getNewmark_Disp_Velocity(ü_n1, u_n, u̇_n, ü_n, t_n1-t_n, β1, β2)

    if isa(f, Function)
        f = f(u_n1, u̇_n1, t_n1)
    end

    fTotal = f + M_vecfunc(ü_n1) + C_vecfunc(u̇_n1) + K_vecfunc(u_n1)

    return fTotal, u_n1, u̇_n1
end

"""The most general form of the Generalized-Alpha method. Gives as output u_n+1, u̇_n+1, u_n+1-αf, u̇_n+1-αf, ü_n+1-αm, αf, αm, γ, and β.
u̇_n+1, u_n+1-αf, u̇_n+1-αf, ü_n+1-αm are the velocity, displacement and acceleration at the end of the time step minus the αf factor for 
the displacement and velocity. The acceleration ü_n+1-αm is the acceleration at the end of the time step minus the αm factor. 
The factors αf is the shift for the displacement, velocity and force vectors. αm is the shift for the acceleration vector. 
u_n+1-αf, u̇_n+1-αf, ü_n+1-αm must be used to calculate the mass, damping, and stiffness matrices. u_n+1, u̇_n+1 are for the next time step.
γ and β are the Newmark parameters."""
function getGenAlpha_DispVelAcc(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    t_n1::Real, t_n::Real, α_f::Real, α_m::Real, γ::Real, β::Real)

    u_n1, u̇_n1 = getNewmark_Disp_Velocity(ü_n1, u_n, u̇_n, ü_n, t_n1-t_n, γ, 2.0*β)
    ü_n1_m_αm = α_m*ü_n + (1.0-α_m)*ü_n1
    u̇_n1_m_αf = α_f*u̇_n + (1.0-α_f)*u̇_n1
    u_n1_m_αf = α_f*u_n + (1.0-α_f)*u_n1
    return u_n1, u̇_n1, u_n1_m_αf, u̇_n1_m_αf, ü_n1_m_αm, α_f, α_m, γ, β
end

"""Given the density parameter ρInf, this function returns the implicit Generalized-Alpha parameters αf, αm, γ, and β.
See Explicit time integration algorithms for structural dynamics with optimal numerical dissipation, J. Chung and G. M. Hulbert, 1996."""
function getImplicitGenAlphaParameters(ρInf::Real)
    α_m = (2.0*ρInf - 1.0)/(ρInf + 1.0)
    α_f = ρInf/(ρInf + 1.0)
    γ = 0.5 - α_m + α_f
    β = 0.25*(0.5 + γ)^2
    return α_f, α_m, γ, β
end

"""Given the density parameter ρb, this function returns the explicit Generalized-Alpha parameters αf, αm, γ, and β.
See Explicit time integration algorithms for structural dynamics with optimal numerical dissipation, J. Chung and G. M. Hulbert, 1996."""
function getExplicitGenAlphaParameters(ρb::Real)
    α_m = (2.0*ρb - 1.0)/(ρb + 1.0)
    β = (5.0 - 3.0*ρb)/((1.0 + ρb)^2.0 * (2.0 - ρb))
    γ = 3.0/2.0 - α_m
    return 1.0, α_m, γ, β
end

"""Given the density parameter ρ and the method, this function returns the Generalized-Alpha parameters αf, αm, γ, and β. method can be either
:implicit or :explicit."""
function getGenAlphaParameters(ρ::Real, method::Symbol)
    if method == :implicit
        return getImplicitGenAlphaParameters(ρ)
    elseif method == :explicit
        return getExplicitGenAlphaParameters(ρ)
    else
        error("Invalid method for Generalized-Alpha. Choices are either :implicit or :explicit.")
    end
end

"""Given the density parameter ρ and the method, this function returns the Generalized-Alpha data u_n+1, u̇_n+1, u_n+1-αf, u̇_n+1-αf, ü_n+1-αm, αf, αm, γ, and β.
u̇_n+1, u_n+1-αf, are the velocity and displacement at the end of the time . _n+1-αf, u̇_n+1-αf, are the velocity and displacement at the end of the time step 
minus the αf factor. ü_n+1-αm is the acceleration at the end of the time step minus the αm factor. αf is the shift for the displacement, velocity and force vectors.
αm is the shift for the acceleration vector. u_n+1-αf, u̇_n+1-αf, ü_n+1-αm must be used to calculate the mass, damping, and stiffness matrices. 
u_n+1, u̇_n+1 are for the next time step. γ and β are the Newmark parameters. method can be either :implicit or :explicit."""
function getGenAlpha_DispVelAcc(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    t_n1::Real, t_n::Real, ρ::Real, method::Symbol)

    α_f, α_m,  γ, β = getGenAlphaParameters(ρ, method)
    return getGenAlpha_DispVelAcc(ü_n1, u_n, u̇_n, ü_n, t_n1, t_n, α_f, α_m, γ, β)
end

"""Given the parameter α_f, this function returns the Generalized-Alpha parameters αf, αm, γ, and β."""
function getHhtAlphaParameters(α_f::Real)
    @assert 0.0 <= α_f <= 1.0/3.0 "α_f must be between 0 and 1/3."
    α_m = 1.0
    β = 0.25*(1.0 + α_f)^2
    γ = 0.5 + α_f
    return α_f, α_m, γ, β
end

"""Given the parameter α_f, this function returns the Generalized-Alpha data u_n+1, u̇_n+1, u_n+1-αf, u̇_n+1-αf, ü_n+1-αm, αf, αm, γ, and β."""
function getHhtAlpha_DispVecAcc(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    t_n1::Real, t_n::Real, α_f::Real)

    α_f, α_m, γ, β = getHhtAlphaParameters(α_f)
    return getGenAlpha_DispVelAcc(ü_n1, u_n, u̇_n, ü_n, t_n1, t_n, α_f, α_m, γ, β)
end
    
    
    



""""This function returns the total force acting on the system, given the current acceleration guess in a Newton-Raphson iteration.
    Additionally, it returns the next set of possible displacement, u_n+1, and velocity, u̇_n+1. 
The mass, damping, and stiffness matrices are represented by the functions M_vecfunc(ü_n1), C_vecfunc(u̇_n1), and K_vecfunc(u_n1), respectively. The
external force acting on the system is represented by f. The current displacement, velocity, and acceleration are represented by u_n,
u̇_n, and ü_n, respectively. The time is represented by t_n and t_n1 (t_n+1). The Generalized-Alpha parameters are represented by ρ and method.
method can be either :implicit or :explicit. The external force acting on the system can be a function of the displacement, velocity, and time. 
The function f is the external force and must be defined as f(u, u̇, t)."""
function GeneralizedAlpha(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f::T, u_n::Union{Real, Vector}, 
    u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, t_n1::Real, t_n::Real, ρ::Real, method::Symbol) where T

    u_n1, u̇_n1, u_n1_m_αf,  u̇_n1_m_αf, ü_n1_m_αm, αf, αm, γ, β = getGenAlpha_DispVelAcc(ü_n1, u_n, u̇_n, ü_n, t_n1, t_n, ρ, method)
    
    if isa(f, Function)
        f = αf*f(u_n, u̇_n, t_n)+ (1.0 - αf)*f(u_n1, u̇_n1, t_n1)
    end
    fTotal = f + M_vecfunc(ü_n1_m_αm) + C_vecfunc(u̇_n1_m_αf) + K_vecfunc(u_n1_m_αf)
    return fTotal, u_n1, u̇_n1
end

"""This function returns the total force acting on the system, given the current acceleration guess in a Newton-Raphson iteration
 along with the next set of possible displacement, u_n+1, and velocity, u̇_n+1. The external forces acting on the system are 
 represented by f_n_ext and f_n1_ext at the last and current time steps, respectively. The mass, damping, and stiffness matrices are
represented by the functions M_vecfunc(ü_n1), C_vecfunc(u̇_n1), and K_vecfunc(u_n1), respectively. The current displacement, velocity, and acceleration
are represented by u_n, u̇_n, and ü_n, respectively. The time is represented by t_n and t_n1 (t_n+1). 
The Generalized-Alpha parameters are represented by ρ and method."""
function GeneralizedAlpha(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f_n_ext::Union{Real, Vector},
    f_n1_ext::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, t_n1::Real, t_n::Real, ρ::Real, 
    method::Symbol)

    u_n1, u̇_n1, u_n1_m_αf,  u̇_n1_m_αf, ü_n1_m_αm, αf, αm, γ, β = getGenAlpha_DispVelAcc(ü_n1, u_n, u̇_n, ü_n, t_n1, t_n, ρ, method)
    
    f = αf*f_n_ext + (1.0 - αf)*f_n1_ext

    fTotal = f + M_vecfunc(ü_n1_m_αm) + C_vecfunc(u̇_n1_m_αf) + K_vecfunc(u_n1_m_αf)
    return fTotal, u_n1, u̇_n1
end

