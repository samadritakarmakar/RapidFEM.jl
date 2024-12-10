""" Time integration test functions for 2nd order ODEs. t is the time, m is the mass, c is the damping coefficient, 
k is the spring constant, f is the external force. D, E, and G are the chosen polynomial coefficients, 
The analytical solution is x = D*t^2 + E*t + G. """

 function analyticalForceFunction(t::Float64, m::Float64, c::Float64, k::Float64, D::Float64, E::Float64, G::Float64)
    return m*2D + c*(2*D*t + E) + k*(D*t^2 + E*t + G) #k*D*t^2 + (2*D*c + E*k)*t + (2*D*m + E*c + G*k)
end