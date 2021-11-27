#====================================================================
  Copyright (c) 2020 Samadrita Karmakar samadritakarmakar@gmail.com

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 =====================================================================#
 
include("gaussLine.jl")
function gaussQuad(order::Int64)
    w_base::Array{Float64,1}, ip_base::Array{Float64,1} = gaussLine(order)
    ip_length::Int64 = length(w_base)
    w::Array{Float64,1} = Array{Float64,1}(undef, ip_length^2)
    ip::Array{Float64,2} = Array{Float64,2}(undef, ip_length^2, 2)
    ip_no::Int64 = 1
    for i ∈ 1:ip_length
        for j ∈ 1:ip_length
            w[ip_no] = w_base[i]*w_base[j]
            ip[ip_no, :] = [ip_base[i] ip_base[j]]
            ip_no += 1
        end
    end
    return w, ip
end
