module DroughtIndices()

# Write your package code here.
include("moms.jl")
include("thornthwaite.jl")
include("cdfglo.jl")
include("params.jl")

function SPEI_helper(x_mon, nmom)
    lmoms = moms(x_mon, nmom)
    lambdas = lmoms[1]
    ratios = lmoms[2]
    para = params(lambdas, ratios)
    spei = cdfglo(x_mon, para) 
    return(spei)
end

function calculate_SPEI(Precip, Tave, lat, start_month = 1)
    nmom = 3
    v_dims = size(Precip)[1]
    i_p = deepcopy(Precip)
    Bal = Precip - thornthwaite(Tave, lat, start_month)
    rep_len = (Int64.(v_dims) / 12)
    rep_len = floor.(Int, rep_len)
    remain_len = v_dims % 12
    months = repeat(1:12, rep_len)
    append!(months, months[1:remain_len])
    mat = hcat(months, Bal)
    spei = deepcopy(mat[:,1])
    for i in 1:12
        x_mon = deepcopy(mat[:,2][findall(mat[:,1] .== i)])
        h = SPEI_helper(x_mon, nmom)
        spei[findall(mat[:,1] .== i)] .= (h)
    end
    return(spei)
end

function calculate_SPI(Precip)
    nmom = 3
    v_dims = size(Precip)[1]
    rep_len = (Int64.(v_dims) / 12)
    rep_len = floor.(Int, rep_len)
    remain_len = v_dims % 12
    months = repeat(1:12, rep_len)
    append!(months, months[1:remain_len])
    mat = hcat(months, Precip)
    spi = deepcopy(mat[:,1])
    for i in 1:12
        x_mon = deepcopy(mat[:,2][findall(mat[:,1] .== i)])
        h = SPEI_helper(x_mon, nmom)
        spi[findall(mat[:,1] .== i)] .= (h)
    end
    return(spi)
end
export calculate_SPEI
export calculate_SPI
export thornthwaite
export SPEI_helper
export moms
export params
export cdfglo
end