############################
## Extended Kalman Filter ##
############################

## Model structure
make.model = function(y, n = 1, k = NULL, p = 0, V, W, Wx, omega = 2*pi/365.25,
                      dn = FALSE, Bn = NULL, Zn = NULL, dp = FALSE, Bp = NULL) {

    ## Dimensions
    tt = length(y)                           ## Length of observation vector
    kk = length(k)                           ## Number of seasonal harmonics
    zz = if (is.null(Zn)) 0 else ncol(Zn)    ## Number of covariates
    nn = n + 2*kk + 2*p + (dn != 0) + (p+1)*(dp != 0) + zz

    ## Model structure
    rho = 0
    na  = rho + 1; no = na + n    - 1; rho = no
    ka  = rho + 1; ko = ka + 2*kk - 1; rho = ko
    xa  = rho + 1; xo = xa + p    - 1; rho = xo
    pa  = rho + 1; po = pa + p    - 1; rho = po
    if (dn != 0) {
        dna = rho + 1; dno = dna; rho = dno
    } else {
        dna = NULL   ; dno = NULL
    }
    if (!is.null(Zn)) {
        zna = rho + 1; zno = zna + zz - 1; rho = zno
    } else {
        zna = NULL; zno = NULL
    }
    if (dp != 0) {
        dpa = rho + 1; dpo = dpa + p; rho = dpo
    } else {
        dpa = NULL   ; dpo = NULL
    }

    ## Initialise storage
    F  = array(0, dim = c(nn, 1))
    G  = array(0, dim = c(nn,nn))
    Fj = array(0, dim = c( 1,nn))
    Gj = array(0, dim = c(nn,nn))
    Vj = array(0, dim = c( 1, 1))
    Wj = array(0, dim = c(nn,nn))

    ## Polynomial trend component
    if (n > 0) {

        F[na,1] = 1

        J = 1*upper.tri(matrix(0, n, n), diag = TRUE)
        G[na:no,na:no] = J

    }

    ## Seasonal component
    if (kk > 0) {

        F[seq(ka,ko,2),1] = 1

        J = matrix(0, 2*kk, 2*kk)
        for (i in 1:kk)
            J[(2*i-1):(2*i),(2*i-1):(2*i)] =
            c( cos(k[i]*omega),-sin(k[i]*omega),
               sin(k[i]*omega), cos(k[i]*omega))
        G[ka:ko,ka:ko] = J

    }

    ## Latent auto-regressive component
    if (p > 0) {

        F[xa,1] = 1

        if (p > 1) {
            J = diag(p-1)
            G[(xa+1):xo,xa:(xo-1)] = J
        }

        J = diag(p)
        G[pa:po,pa:po] = J

    }

    ## Trend intervention
    if (dn != 0) {

        J = dn
        G[dna:dno,dna:dno] = J

    }

    ## Trend covariates
    if (!is.null(Zn)) {

        J = diag(zz)
        G[zna:zno,zna:zno] = J

    }

    ## Auto-regressive intervention
    if (dp != 0) {

        J = 1
        G[dpa,xo] = 1

        J = diag(rep(dp, p))
        G[(dpa+1):dpo,(dpa+1):dpo] = J

    }

    ## Jacobians
    Fj    = t(F)
    Gj    = G
    Vj[,] = 1
    Wj[,] = 1 * (W > 0)

    Wj[na:no,na:no] = G[na:no,na:no]

    ## Make observation function
    make.Fn = function(F, Fj, Bn, dna, Zn, zna, zno, Bp, dpa, dpo, xa, xo) {
        dn = ! is.null(dna)
        zn = ! is.null(Zn)
        dp = ! is.null(dpa)
        if (zn)
            zna.zno = zna:zno
        if (dp) {
            dpa.dpo = (dpa+1):(dpo-1)
             xa.xo  = (xa +1): xo
        }
        function(a, v, t) {

            ## Modify observation matrix and Jacobians
            if (dn)
                F[dna,1] = Fj[1,dna] = Bn(t)
            if (zn)
                F[zna.zno,1] = Fj[1,zna.zno] = Zn[t,]
            if (dp) {
                bp = Bp(t)
                F[xa.xo,1] = Fj[1,xa.xo] = bp * a[dpa.dpo]
                F[dpa  ,1] = Fj[1,dpa  ] = bp * a[dpo    ]
                Fj[1,dpa.dpo] = bp * a[xa.xo]
                Fj[1,dpo    ] = bp * a[dpa  ]
            }

            ## Compute forecast
            f = crossprod(F, a)

            ## Return results
            return(list(f = f, Fj = Fj))

        }

    }
    Fn = make.Fn(F, Fj, Bn, dna, Zn, zna, zno, Bp, dpa, dpo, xa, xo)
    Fn = cmpfun(Fn)


    ## Make observation variance function
    make.Vn = function(V, Vj)
        function(a, v, t) {

            ## Return results
            return(list(V = V, Vj = Vj))

        }
    Vn = make.Vn(V, Vj)
    Vn = cmpfun(Vn)


    ## Make state function
    make.Gn = function(G, Gj, xa, xo, pa, po) {

        ## Sequences
        xa.xo = xa:xo
        pa.po = pa:po

        function(m, w, t) {

            ## Basic forecast
            a = G %*% m

            ## X
            a [xa      ] = a[pa.po] %*% m[xa.xo]
            Gj[xa,xa.xo] = a[pa.po]
            Gj[xa,pa.po] = m[xa.xo]

            ## Return result
            return(list(a = a, Gj = Gj))

        }

    }
    Gn = make.Gn(G, Gj, xa, xo, pa, po)
    Gn = cmpfun(Gn)


    ## Make state variance function
    make.Wn = function(W, Wj, Wx, xa, xo, pa, po) {

        ## Sequences
        xa.xo   = xa:xo
        pa.po   = pa:po

        function(m, w, t) {

            ## Modify state variance
            W[xa,xa] = Wx(t)

            ## X
            Wj[xa,pa.po] = m[xa.xo]

            ## Return result
            return(list(W = W, Wj = Wj))

        }

    }
    Wn = make.Wn(W, Wj, Wx, xa, xo, pa, po)
    Wn = cmpfun(Wn)

    ## Return results
    return(list(Fn = Fn, Gn = Gn, Vn = Vn, Wn = Wn))

}
