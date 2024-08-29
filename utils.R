
## given CDF, density and quantile function, return same after left truncation
trunc.distr <- function(lower,p.full,d.full,q.full) {
    p.trunc <- function(q) (q>=lower)*(p.full(q)-p.full(lower))/(1-p.full(lower))
    d.trunc <- function(x) (x>=lower)*d.full(x)/(1-p.full(lower))
    q.trunc <- function(p) q.full(p*(1-p.full(lower)) + p.full(lower))
    r.trunc <- function(n) q.trunc(runif(n))
    return(list(p=p.trunc,d=d.trunc,q=q.trunc,r=r.trunc))
}


Z.unif <- function(supp.Z=c(-sqrt(12)/2,sqrt(12)/2)) {
    cutoff <- supp.Z[1]
    E.f <- 1/(supp.Z[2]*2)
    E.ZF <-  1/(supp.Z[2]*2)^2 * 1/3*((supp.Z[2]*2)^3/4) #+ s/4*(s^2/4-cutoff^2))
    var.Z <- 1/12 * (supp.Z[2]*2)^2
    rZ <- function(n) runif(n,cutoff,supp.Z[2])
    dZ <- function(z) dunif(z,cutoff,supp.Z[2])
    pZ <- function(q)punif(q,cutoff,supp.Z[2])
    theta.to.cutoff <- function(theta)2*theta-supp.Z[2]
    return( structure(class='Z.distribution',mget(ls())) )
}

Z.normal <- function(supp.Z=c(-Inf,Inf)) {
    cutoff <- supp.Z[1]
    normal.trunc <- trunc.distr(lower=cutoff,p.full=pnorm,d.full=dnorm,q.full=qnorm)
    rZ <- normal.trunc$r
    dZ <- normal.trunc$d
    pZ <- normal.trunc$p
    qZ <- normal.trunc$q
    theta.to.cutoff <- function(theta) {
        if(theta==0)return(-Inf)
        with(list(    mu = function(x)exp(log(dnorm(x))-log(1-pnorm(x))) ),
             uniroot(function(x)mu(x)-theta,c(-1,1),extendInt='yes')$root  )
    }
    E.f <- E.ZF <- 1/(2*sqrt(pi))
    var.Z <- 1
    return( structure(class='Z.distribution',mget(ls())) )
}


Z.t <- function(df) {
    stopifnot(df>2)
    function(supp.Z=c(-Inf,Inf)) {
        cutoff <- supp.Z[1]
        sigma <- sqrt(df/(df-2))
        const <- gamma((df+1)/2)/sqrt(pi*df)/gamma(df/2) ## integrating factor
        t.trunc <- trunc.distr(lower=cutoff,p.full=function(q)pt(q*sigma,df),
                               d.full=function(x)dt(x*sigma,df)*sigma,
                               q.full=function(p)1/sigma*qt(p,df))
        rZ <- t.trunc$r
        dZ <- t.trunc$d
        pZ <- t.trunc$p
        qZ <- t.trunc$q
        theta.to.cutoff <- function(theta) {
            if(theta==0)return(-Inf)
            mu  <-  function(cutoff)1/(1-pt(cutoff*sigma,df))*const/sigma*df/(df-1)*(1+sigma^2*cutoff^2/df)^((1-df)/2)
            uniroot(function(x)mu(x)-theta,c(-1,1),extendInt='yes')$root  
        }
        E.f <-1/sqrt(pi)/sqrt(df-2)*gamma((df+1)/2)^2/gamma(df/2)^2*gamma(df+1/2)/gamma(df+1)
        E.ZF <- 2*const^2*df^2*sqrt(df-2)/(df-1)*beta(3/2,df-1/2)
        var.Z <- 1
        return( structure(class='Z.distribution',mget(ls())) )
    }
}


Z.beta <- function(a,b) {
    stopifnot(a>0 && b>0)
    function(supp.Z= (c(0,1) - a/(a+b)) / sqrt(a*b/(a+b)^2/(a+b+1))) {
        cutoff <- supp.Z[1]
        mu <- a/(a+b)
        sigma2 <- a*b/(a+b)^2/(a+b+1)
        sigma <- sqrt(sigma2)
        beta.trunc <- trunc.distr(lower=cutoff,p.full=function(q)pbeta(q*sigma+mu,a,b),
                               d.full=function(x)dbeta(x*sigma+mu,a,b)*sigma,
                               q.full=function(p)(qbeta(p,a,b)-mu)/sigma)
        rZ <- beta.trunc$r
        dZ <- beta.trunc$d
        pZ <- beta.trunc$p
        qZ <- beta.trunc$q
        theta.to.cutoff <- function(theta) {
            if(theta==0)return(-Inf)
            cutoff.to.theta  <-  function(cutoff)-mu/sigma + 1/sigma/(1-pbeta(sigma*cutoff+mu,a,b))*beta(a+1,b)/beta(a,b)*(1-pbeta(cutoff*sigma+mu,a+1,b))
            uniroot(function(x)cutoff.to.theta(x)-theta,c(-1,1),extendInt='yes')$root  
        }
        E.f <- sigma*beta(2*a-1,2*b-1)/beta(a,b)^2
        E.ZF <- with(list(d.full=function(x)dbeta(x*sigma+mu,a,b)*sigma,
                          p.full=function(q)pbeta(q*sigma+mu,a,b),
                          supp.Z= (c(0,1) - a/(a+b)) / sqrt(a*b/(a+b)^2/(a+b+1))),
                     integrate(function(z)z*d.full(z)*p.full(z),supp.Z[1],supp.Z[2])$val)
        var.Z <- 1
        return( structure(class='Z.distribution',mget(ls())) )
    }
}



S.unif <- function(S.par=c(0,1)) {
    supp.S <- S.par
    dS <- function(s)dunif(s,supp.S[1],supp.S[2])
    pS <- function(q)punif(q,supp.S[1],supp.S[2])
    rS <- function(n)runif(n,supp.S[1],supp.S[2])
    E.S2 <- diff(supp.S)^2/12 + mean(supp.S)^2
    E.S1 <- mean(supp.S)
    mu.S <- function(k)integrate(function(s)s^k*dS(s),supp.S[1],supp.S[2])$value
    mean.S.pair <- with(list(b=supp.S[2],a=supp.S[1]), 1/(b-a)*((a^2+b^2)/3-2/3*a*b))
    var.S <- 1/12*diff(supp.S)^2
    return( structure(class='S.distribution',mget(ls())) )
}


S.beta <- function(S.par=c(a=1,b=1,shift=0)) {
    dS <- function(s)dbeta(s-S.par['shift'],S.par['a'],S.par['b'])
    pS <- function(q)pbeta(q-S.par['shift'],S.par['a'],S.par['b'])
    rS <- function(n)rbeta(n,S.par['a'],S.par['b'])+S.par['shift']
    E.S1 <- with(as.list(S.par),a/(a+b)) + S.par['shift']
    var.S <- with(as.list(S.par), a*b/(a+b)^2/(a+b+1))
    supp.S <- c(0,1)+S.par['shift']
    E.S2 <- E.S1^2 + var.S
    mean.S.pair <- 2 * integrate(function(s)pS(s)*(1-pS(s)), supp.S[1],supp.S[2])$val
    return( structure(class='S.distribution',mget(ls())) )
}


S.lognormal <- function(S.par=c(meanlog=0,sdlog=1,shift=0)) {
    dS <- function(s)dlnorm(s-S.par['shift'],S.par['meanlog'],S.par['sdlog'])
    pS <- function(q)plnorm(q-S.par['shift'],S.par['meanlog'],S.par['sdlog'])
    rS <- function(n)rlnorm(n,S.par['meanlog'],S.par['sdlog'])+S.par['shift']
    E.S1 <- with(as.list(S.par),exp(meanlog+sdlog^2/2)) + S.par['shift']
    var.S <- with(as.list(S.par),  (exp(sdlog^2)-1)*exp(2*meanlog+sdlog^2) )
    supp.S <- c(0,Inf)+S.par['shift']
    E.S2 <- E.S1^2 + var.S
    mean.S.pair <- 2 * integrate(function(s)pS(s)*(1-pS(s)), supp.S[1],supp.S[2])$val
    return( structure(class='S.distribution',mget(ls())) )
}

## type 1 pareto. EnvStats::pareto man page maps to wikipedia
## paremters as shape=theta --> alpha, location=eta --> sigma
S.pareto <- function(S.par=c(location=1,shape=1)) {
    stopifnot(S.par['shape'] > 2)
    dS <- function(s)EnvStats::dpareto(s,location=S.par['location'],shape=S.par['shape'])
    pS <- function(q)EnvStats::ppareto(q,location=S.par['location'],shape=S.par['shape'])
    rS <- function(s)EnvStats::rpareto(n,location=S.par['location'],shape=S.par['shape'])
    E.S1 <- with(as.list(S.par), location*shape / (shape - 1))
    E.S2 <- with(as.list(S.par), location^2*shape / (shape - 2))
    var.S <- E.S2 - (E.S1)^2
    supp.S <- c(0,Inf)+S.par['location']
    mean.S.pair <- 2 * integrate(function(s)pS(s)*(1-pS(s)), supp.S[1],supp.S[2])$val
    return( structure(class='S.distribution',mget(ls())) )
    }


## ## TODO: conform to htest class, also for Egger's
## begg.test <- function(y,v,method='kendall',exact=FALSE,...) {
##     theta.fe <- sum(y/v)/sum(1/v)
##     if(!exact) {
##         out <- with(cor.test((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method=method,...),
##              structure(unname(c(stat=estimate,pval=p.value)),names=c('statistic','p.value')))
##     } else {
##         s <- 1/sqrt(v)
##         z <- y*s
##         theta.fe <- sum(z*s)/sum(s^2) # already calculated, clean
##         n <- length(z)
##         tau.hat <- 2*mean(apply(combn(n,2),2,function(idx)(z[idx[1]]-z[idx[2]])/(s[idx[1]]-s[idx[2]])<theta.fe)) - 1
##         pval <- 2*(1-pnorm(sqrt(9*n/4)*abs(tau.hat)))
##         c(stat=tau.hat,pval=pval)
##     }
## }

begg.test <- function(y,v,method='kendall',exact=FALSE,...) {
    theta.fe <- sum(y/v)/sum(1/v)
    if(!exact) {
        out <- cor.test((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method=method,...)
             ## structure(unname(c(stat=estimate,pval=p.value)),names=c('statistic','p.value')))
    } else {
        s <- 1/sqrt(v)
        z <- y*s
        ## theta.fe <- sum(z*s)/sum(s^2) # already calculated, clean
        n <- length(z)
        tau.hat <- 2*mean(apply(combn(n,2),2,function(idx)(z[idx[1]]-z[idx[2]])/(s[idx[1]]-s[idx[2]])<theta.fe)) - 1
        statistic <- sqrt(9*n/4)*abs(tau.hat)
        pval <- 2*(1-pnorm(statistic))
        out <- list(estimate=tau.hat, p.value=pval, statistic=statistic, stderr = 2/3/sqrt(n), null.value=0, parameter=c(), alternative='two-sided')
    }
    out$data.name <-  paste(deparse(substitute(y)), 'and', deparse(substitute(v)))
    out$method <- paste("Begg's test", if(exact) "(exact)" else "(cor.test approximation)")
    class(out) <- 'htest'
    return(out)
}


## egger.test <- function(y,v,robust=FALSE) {
##     lm0 <- if(robust) {
##                estimatr::lm_robust(I(y/sqrt(v)) ~ I(1/sqrt(v)))
##            } else {
##                lm(I(y/sqrt(v)) ~ I(1/sqrt(v)))
##            }
##     structure(unname(coef(summary(lm0))[1,c(1,3,4)]), names=c('statistic','z.statistic','p.value'))
## }

egger.test <- function(y,v,robust=FALSE) {
    lm0 <- if(robust) {
               estimatr::lm_robust(I(y/sqrt(v)) ~ I(1/sqrt(v)))
           } else {
               lm(I(y/sqrt(v)) ~ I(1/sqrt(v)))
           }
    out <- structure(unname(coef(summary(lm0))[1,]), names=c('estimate','stderr','statistic','p.value'))
    out <- c(out, list(null.value = 0, parameter = c(), alternative = 'two-sided', method = paste0("Egger's test",if(robust)' (robust)' else NULL), data.name = paste(deparse(substitute(y)), 'and', deparse(substitute(v))) ))
    class(out) <- 'htest'
    return(out)
}

ZS.distr <- function(Z.distr,par.S,supp.Z=NULL,cutoff.mean=0) {
    supp.Z[1] <- Z.distr(supp.Z)$theta.to.cutoff(cutoff.mean)
    par.Z <- Z.distr(supp.Z)
    asy.var.begg <- with(par.Z, with(par.S,
                                     4/9 + 4*mean.S.pair^2/E.S2*E.f*(E.f*var.Z-2*E.ZF) 
                         ))
    asy.var.egger <- with(par.Z, with(par.S,
                                      E.S2*var.Z/var.S
                                      ))
    slope.begg <- with(par.Z, with(par.S,
                                     2*E.S1/E.S2*E.f*mean.S.pair / sqrt(asy.var.begg) 
                         ))
    slope.egger <- 1/sqrt(asy.var.egger)
    ARE <- (slope.egger/slope.begg)^2
    structure(class='ZS.distribution', c(par.Z, par.S, mget(ls())))
}

