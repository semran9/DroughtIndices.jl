# using Statistics
function thornthwaite(Tave, lat, start_month = 1)
    v_dims = size(Tave)[1]
    mlen = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    m_count = 1:12
    m_day_helper = cumsum(mlen) .- mlen .+ 15
    m_day = Vector{Any}(undef, v_dims)
    rep_len = (Int64.(v_dims) / 12)
    rep_len = floor.(Int, rep_len)
    remain_len = v_dims % 12
    m_day = cumsum(mlen) .- mlen .+ 15
    m_day = repeat(m_day, rep_len)
    mlen_full = repeat(mlen, rep_len)
    m_count_full = repeat(m_count, rep_len)
    append!(m_day, m_day_helper[1:remain_len])
    append!(mlen_full, mlen[1:remain_len])
    append!(m_count_full, m_count[1:remain_len])
    if start_month != 1 
        end_vec_mlen = mlen_full[1:(start_month-1)] 
        end_vec_mcount = m_count_full[1:(start_month-1)]
        m_count_full = m_count_full[start_month:size(m_count_full)[1]]
        mlen_full = mlen_full[start_month:size(mlen_full)[1]]
        append!(mlen_full, end_vec_mlen)
        append!(m_count_full, end_vec_mcount)
    end
    T_i = deepcopy(Tave[1:(v_dims)])
    lat_i = Vector{Float64}(undef, v_dims)
    lat_i = fill!(lat_i, lat)   
    # correction factor
    tanLat = tan.(lat_i ./ 57.2957795)
    # solar declination
    Delta = 0.4093 .* sin.(((2 * pi .* m_day) ./ 365) .- 1.405)
    # sun rising angle
    tanDelta = tan.(Delta)
    tanLatDelta = (tanLat .* tanDelta)
    tanLatDelta[findall(tanLatDelta .< -1)] .= -1
    tanLatDelta[findall(tanLatDelta .> 1)] .= 1
    omega = acos.(-tanLatDelta)
    # mean daily daylight hours for each month (N)
    N = (24 / pi) .* omega
    # which leads to K
    K = (N / 12) .* (mlen_full ./ 30)
    Tt = Vector{Any}(undef, v_dims)
    for i in 1:min(12, v_dims)
        Tt[i] = mean(T_i[findall(m_count_full .== i)])
    end
    Tt = convert(Vector{Float64}, Tt)
    Tt[findall(Tt .< Float64(0.0))] .= 0.0
    J = sum((Tt ./ 5).^1.514)
    J2 = J*J
    J3 = J2*J
    q =  0.000000675 * J3 - 0.0000771 * J2 + 0.01792 * J + 0.49239
    T_i[findall(T_i .< 0)] .= 0
    if J != 0
        PET = K .* 16 .* (10 .* T_i ./ J).^q
    else
        PET = K .* 0
    end
    return PET
end
export thornthwaite