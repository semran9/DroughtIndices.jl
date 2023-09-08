using  Distributions
function cdfglo(x_mon, para)
    SMALL = 1e-15

    XI = para[1]
    A = para[2]
    K = para[3]

    f = Vector{Float64}(undef, length(x_mon))
    Y = (x_mon.-XI)./A
    if K == 0
        f <- 1 ./ (1 .+ exp.(-Y))
    else 
        ARG = 1 .- K .* Y
        Y = -log.(ARG) ./K
        f = 1 ./ (1 .+ exp.(-Y))
    end
    if K<0
        f[findall(!isfinite, f)] .= 0
    elseif K > 0
        f[findall(!isfinite, f)] .= 1
    end
    return(Distributions.quantile.(Normal(0, 1), f))
end
export cdfglo

