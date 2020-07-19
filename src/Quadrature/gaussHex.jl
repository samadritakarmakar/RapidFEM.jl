include("gaussQuad.jl")
function getQuadratureHex(order::Int64)
    w_base::Array{Float64,1}, ip_base::Array{Float64,1} = getQuadratureLine(order)
    ip_length::Int64 = length(w_base)
    w::Array{Float64,1} = Array{Float64,1}(undef, ip_length^3)
    ip::Array{Float64,2} = Array{Float64,2}(undef, ip_length^3, 3)
    ip_no::Int64 = 1
    for i ∈ 1:ip_length
        for j ∈ 1:ip_length
            for k ∈ 1:ip_length
                w[ip_no] = w_base[i]*w_base[j]*w_base[k]
                ip[ip_no, :] = [ip_base[i] ip_base[j] ip_base[k]]
                ip_no += 1
            end
        end
    end
    return w, ip
end
