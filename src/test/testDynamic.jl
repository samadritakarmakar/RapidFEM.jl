using RapidFEM, SparseArrays

M = spzeros(1,1)
M[1,1] = 1.0
C = spzeros(1,1)
C[1,1] = 2.0
K = spzeros(1,1)
K[1,1] = 3.0
u = zeros(1)
u[1] = 1.0
u̇ = zeros(1)
u̇[1] = 1.0
SolutionArray = [u]
f_Array = [[1.0], [2.0]]
Δt = 0.01
θ_Array = [0.5]
for i ∈ 1:5
A_global, f_mean_global = SSpj_getFinal_A_b([C, K], SolutionArray, f_Array, Δt, θ_Array)
α_global = A_global\f_mean_global
println(" A_global = " , A_global[1], " f_mean_global = ", f_mean_global[1], " α = ", α_global[1])
updateSolution!(SolutionArray, Δt, α_global)
println(SolutionArray, " ")
end
