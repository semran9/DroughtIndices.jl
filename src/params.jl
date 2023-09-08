using Distributions
function params(lambdas, ratios)
    SMALL = 1e-6
    para = Vector{Float64}(undef,3)


    K = -1*ratios[3]
   if abs(K) <= SMALL
        para[3] <- 0
        para[2] <- lambdas[2]
        para[1] <- lambdas[1]
    end
    KK = K*pi/sin(K*pi)
    A  = lambdas[2]/KK
    para[1] = lambdas[1] - A*(1-KK)/K
    para[2] = A
    para[3] = K
    return(para)
end
export params
