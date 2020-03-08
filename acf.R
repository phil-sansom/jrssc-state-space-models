myacf = function(x, lag.max, index) {

    indices = unique(index)
    n.indices = length(indices)

    z = array(NA, c(lag.max+1,n.indices))
    z[1,] = 1

    for (i in 1:n.indices) {

        mask = which(index == indices[i])

        for (l in 1:lag.max) {

            maskl = mask[mask > l]
            z[l+1,i] = cor(x[maskl], x[maskl-l], "complete.obs")

        } ## j

    } ## i

    return(list(lag = 0:lag.max, index = indices, acf = z))

}
