
ar2pacf = function(phi) {
    p = length(phi)
    a = list()
    a[[p]] = phi
    rho = numeric(p)
    rho[p] = a[[p]][p]
    for (j in (p-1):1) {
        for (i in 1:j)
            a[[j]][i] = (a[[j+1]][i] + rho[j+1]*a[[j+1]][j-i+1])/(1-rho[j+1]^2)
        rho[j] = a[[j]][j]
    }
    return(phi)
}

pacf2ar = function(rho) {
    p = length(rho)
    a = list()
    a[[1]] = rho[1]
    for (j in 2:p) {
        a[[j]] = numeric(j)
        for (i in 1:(j-1))
            a[[j]][i] = a[[j-1]][i] - rho[j]*a[[j-1]][j-i]
        a[[j]][j] = rho[j]
    }
    return(a[[p]])
}
