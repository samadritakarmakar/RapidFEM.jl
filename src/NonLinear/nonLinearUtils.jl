"Bergan Increment can be used when in Solid Mechanics problems. If the problem
is force driven and there is a chance of softening behaviour the Bergan factor
can be used to scale the initially applied force vector."
function berganIncrement(f₀::Array{Float64, 1}, P::Array{Float64, 1})
    Δλ = -f₀'*P/(f₀'*f₀)
    return Δλ[1]
end
