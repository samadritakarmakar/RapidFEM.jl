

function polynomial2ndOrder_usol(t::Float64, paramsPolynomial::Tuple)
    m, c, k, D, E, G = paramsPolynomial
    return D*t^2 + E*t + G
end

function polynomial2ndOrder_u̇sol(t::Float64, paramsPolynomial::Tuple)
    m, c, k, D, E, G = paramsPolynomial
    return 2*D*t + E
end

function polynomial2ndOrder_üsol(t::Float64, paramsPolynomial::Tuple)
    m, c, k, D, E, G = paramsPolynomial
    return 2*D
end

""" Time integration test functions for 2nd order ODEs. t is the time, m is the mass, c is the damping coefficient, 
k is the spring constant, f is the external force. D, E, and G are the chosen polynomial coefficients, 
The analytical solution is x = D*t^2 + E*t + G. """
function polynomial2ndOrderF_Func(t::Float64, paramsPolynomial::Tuple)
    m, c, k, D, E, G = paramsPolynomial
    #return m*2D + c*(2*D*t + E) + k*(D*t^2 + E*t + G) #k*D*t^2 + (2*D*c + E*k)*t + (2*D*m + E*c + G*k)
    return m*polynomial2ndOrder_üsol(t, paramsPolynomial) + c*polynomial2ndOrder_u̇sol(t, paramsPolynomial) + k*polynomial2ndOrder_usol(t, paramsPolynomial)
end

polynomialModel2ndOrder = (polynomial2ndOrder_usol, polynomial2ndOrder_u̇sol, polynomial2ndOrder_üsol, polynomial2ndOrderF_Func) 


function oscillator_usol(t::Float64, params::Tuple)
    m, c, k, D = params
    return sin(2*π*t/D) + t^2
end

function oscillator_u̇sol(t::Float64, params::Tuple)
    m, c, k, D  = params
    return 2*π/D*cos(2*π*t/D) + 2*t
end

function oscillator_üsol(t::Float64, params::Tuple)
    m, c, k, D  = params
    return -(2*π/D)^2*sin(2*π*t/D) + 2
end

function oscillator_F_Func(t::Float64, params::Tuple)
    m, c, k, D = params
    return m*oscillator_üsol(t, params) + c*oscillator_u̇sol(t, params) + k*oscillator_usol(t, params)
end

oscillatorModel = (oscillator_usol, oscillator_u̇sol, oscillator_üsol, oscillator_F_Func)

function rampFlat_üsol(t::Float64, params::Tuple)
    m, c, k, u_total, TaccEnd, Ttotal = params

    u_velFinal = u_total/(0.5*TaccEnd + (Ttotal - TaccEnd))

    if t > TaccEnd
        return 0.0
    end
    return u_velFinal/TaccEnd
end


function rampFlat_u̇sol(t::Float64, params::Tuple)
    m, c, k, u_total, TaccEnd, Ttotal = params

    u_velFinal = u_total/(0.5*TaccEnd + (Ttotal - TaccEnd))

    if t > TaccEnd
        return u_velFinal
    end
    return u_velFinal*t/TaccEnd
end

function rampFlat_usol(t::Float64, params::Tuple)
    m, c, k, u_total, TaccEnd, Ttotal = params

    u_velFinal = u_total/(0.5*TaccEnd + (Ttotal - TaccEnd))
    if t > TaccEnd
        return 0.5*u_velFinal*TaccEnd + u_velFinal*(t - TaccEnd)
    end
    return 0.5*u_velFinal*t^2/TaccEnd
end

function rampFlat_F_Func(t::Float64, params::Tuple)
    m, c, k, u_total, TaccEnd, Ttotal = params
    return m*rampFlat_üsol(t, params) + c*rampFlat_u̇sol(t, params) + k*rampFlat_usol(t, params)
end

rampFlatModel = (rampFlat_usol, rampFlat_u̇sol, rampFlat_üsol, rampFlat_F_Func)




