function moms(x_mon, nmom)
    nmom = 3
    betas = Vector{Float64}(undef, nmom) * 0
    n = length(x_mon)
    x_mon1 = sort(x_mon)
    for r in 0:(nmom-1)
        sum_t = 0.0
        i = r+1
        for j in 1:n
            sum_t +=  binomial(j-1, r) * x_mon1[j] 
        end
        betas[i] = sum_t / (n*binomial(n-1, r))
    end
    lambdas = Vector{Float64}(undef, nmom)
    ratios = deepcopy(lambdas)
    for i in 1:nmom
        r = i - 1
        sum = 0
        for k in 0:r 
            weight = (-1)^(r-k)*binomial(r,k)*binomial(r+k,k)
            sum = sum + weight*betas[k+1]
            lambdas[i] = sum
        end
        if nmom >= 2
        ratios[2] = lambdas[2]/lambdas[1]
        end
        if nmom >= 3
        for r in 3:nmom
            ratios[r] = lambdas[r]/lambdas[2]
        end
        end
    end
    output = []
    push!(output, lambdas)
    push!(output, ratios)
    return(output)
end
export moms