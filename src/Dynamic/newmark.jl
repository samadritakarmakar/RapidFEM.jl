"""For non-linear problems, this function returns the next set of possible displacement, velocity, given the current acceleration guess in 
a Newton-Raphson iteration. WARNING: β2 parameter is here is 2*β parameter in the literature. β1 = γ as in the literature."""
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
function ImplicitNewmark(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f::T, u_n::Union{Real, Vector}, 
    u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, t_n1::Real, t_n::Real, β1::Real, β2::Real) where T

    u_n1, u̇_n1 = getNewmark_Disp_Velocity(ü_n1, u_n, u̇_n, ü_n, t_n1-t_n, β1, β2)

    if isa(f, Function)
        f = f(u_n1, u̇_n1, t_n1)
    end

    fTotal = f + M_vecfunc(ü_n1) + C_vecfunc(u̇_n1) + K_vecfunc(u_n1)

    return fTotal, u_n1, u̇_n1
end

function getGenAlpha_DispVelAccTime(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    t_n1::Real, t_n::Real, α_f::Real, α_m::Real, γ::Real, β::Real)

    u_n1, u̇_n1 = getNewmark_Disp_Velocity(ü_n1, u_n, u̇_n, ü_n, t_n1-t_n, γ, 2.0*β)
    ü_n1_m_αm = α_m*ü_n + (1.0-α_m)*ü_n1
    u̇_n1_m_αf = α_f*u̇_n + (1.0-α_f)*u̇_n1
    u_n1_m_αf = α_f*u_n + (1.0-α_f)*u_n1
    return u_n1, u̇_n1, u_n1_m_αf, u̇_n1_m_αf, ü_n1_m_αm, α_f
end

function getImplicitGenAlphaParameters(ρInf::Real)
    α_m = (2.0*ρInf - 1.0)/(ρInf + 1.0)
    α_f = ρInf/(ρInf + 1.0)
    γ = 0.5 - α_m + α_f
    β = 0.25*(0.5 + γ)^2
    return α_f, α_m, γ, β
end

function getExplicitGenAlphaParameters(ρb::Real)
    α_m = (2.0*ρb - 1.0)/(ρb + 1.0)
    β = (5.0 - 3.0*ρb)/((1.0 + ρb)^2.0 * (2.0 - ρb))
    γ = 3.0/2.0 - α_m
    return 1.0, α_m, γ, β
end

function getGenAlphaParameters(ρ::Real, method::Symbol)
    if method == :implicit
        return getImplicitGenAlphaParameters(ρ)
    elseif method == :explicit
        return getExplicitGenAlphaParameters(ρ)
    else
        error("Invalid method for Generalized-Alpha. Choices are either :implicit or :explicit.")
    end
end

function getGenAlpha_DispVelAccTime(ü_n1::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, 
    t_n1::Real, t_n::Real, ρ::Real, method::Symbol)

    α_f, α_m,  γ, β = getGenAlphaParameters(ρ, method)
    return getGenAlpha_DispVelAccTime(ü_n1, u_n, u̇_n, ü_n, t_n1, t_n, α_f, α_m, γ, β)
end

function GeneralizedAlpha(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f::T, u_n::Union{Real, Vector}, 
    u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, t_n1::Real, t_n::Real, ρ::Real, method::Symbol) where T

    u_n1, u̇_n1, u_n1_m_αf,  u̇_n1_m_αf, ü_n1_m_αm, αf = getGenAlpha_DispVelAccTime(ü_n1, u_n, u̇_n, ü_n, t_n1, t_n, ρ, method)
    
    if isa(f, Function)
        f = αf*f(u_n, u̇_n, t_n)+ (1.0 - αf)*f(u_n1, u̇_n1, t_n1)
    end
    fTotal = f + M_vecfunc(ü_n1_m_αm) + C_vecfunc(u̇_n1_m_αf) + K_vecfunc(u_n1_m_αf)
    return fTotal, u_n1, u̇_n1
end

function GeneralizedAlpha(ü_n1::Union{Real, Vector}, M_vecfunc::Function, C_vecfunc::Function, K_vecfunc::Function, f_n_ext::Union{Real, Vector},
    f_n1_ext::Union{Real, Vector}, u_n::Union{Real, Vector}, u̇_n::Union{Real, Vector}, ü_n::Union{Real, Vector}, t_n1::Real, t_n::Real, ρ::Real, 
    method::Symbol)

    u_n1, u̇_n1, u_n1_m_αf,  u̇_n1_m_αf, ü_n1_m_αm, αf = getGenAlpha_DispVelAccTime(ü_n1, u_n, u̇_n, ü_n, t_n1, t_n, ρ, method)
    
    f = αf*f_n_ext + (1.0 - αf)*f_n1_ext

    fTotal = f + M_vecfunc(ü_n1_m_αm) + C_vecfunc(u̇_n1_m_αf) + K_vecfunc(u_n1_m_αf)
    return fTotal, u_n1, u̇_n1
end

