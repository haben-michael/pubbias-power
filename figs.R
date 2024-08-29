
## 1. figure illustrating beggs test


## Y/sigma coordinates
require(MASS)
n <- 1e6
extrafont::loadfonts()
png('figs/thresholding_y.png', width = 1024, height = 768, pointsize=15, family='CM Roman')
threshold <- qnorm(1-.3)
sigma <- sort(runif(n,0,1))
y <- rnorm(n,sd=sigma)
smoothScatter(sigma,y,colramp=function(n)gray.colors(n,rev=TRUE),ylab='Y',xlab=expression(sigma),xlim=c(.1,max(sigma)-.2),ylim=c(-2,2),bandwidth=.1,transformation=function(x)x^.25)
abline(h=0,lty=1)
abline(0,threshold,lty=2)
y.conditional.mean <- dnorm(threshold)/(1-pnorm(threshold))*sigma
z.conditional.mean <- dnorm(threshold)/(1-pnorm(threshold))
abline(a=0,b=z.conditional.mean,lty=3)
dev.off()

## Z/S coordinates
png('figs/thresholding_z.png', width = 1024, height = 768, pointsize=15, family='CM Roman')
z <- y/sigma
smoothScatter(sigma,z,colramp=function(n)gray.colors(n,rev=TRUE),ylab='Z',xlab='S',xlim=c(.1,max(sigma)-.2),ylim=c(-2,2),bandwidth=.1,transformation=function(x)x^.25)
abline(h=0,lty=1)
abline(h=threshold,lty=2)
z.conditional.mean <- dnorm(threshold)/(1-pnorm(threshold))
abline(h=z.conditional.mean,lty=3)
abline(a=z.conditional.mean,b=-z.conditional.mean,lty=4)
dev.off()




## 2. power simulations

source('utils.R')
sim.power <- function(Z.distr,par.S,theta,alpha, B=B,n=n) {
    ZS <- ZS.distr(Z.distr=Z.distr,par.S=par.S,cutoff.mean=theta) 
    p.vals <- replicate(B,  expr={
        s <- ZS$rS(n)
        z <- ZS$rZ(n)
        y <- z/s
        if(any(is.infinite(y)))browser()
        v <- 1/s^2
        p.val.egger <- egger.test(y,v)['p.value']
        p.val.begg <- begg.test(y,v)['p.value']
        c(egger=unname(p.val.egger), begg=unname(p.val.begg))
    })
    rowMeans(p.vals < alpha)
}
sim.plot <- function(by.params,split.par,thetas,filename=NULL) {
    if(!is.null(filename)) {
        print(filename)
        extrafont::loadfonts()
        png(paste0('figs/',filename), width = 1024, height = 768, pointsize=15, family='CM Roman')
    }
    op <- par(mfrow=c(1,3))
    for(test in c('egger','begg')) {
        label <- switch(test,egger="Egger's test",begg="Begg's test")
        plot(0,type='n',ylim=c(0,1),xlim=range(thetas),xlab=expression(theta),ylab='power',main=label)
        abline(h=alpha,lty=2)
        for(df in split(by.params,by.params[,split.par])) {
            lines(df$theta, df[,test])
            text(df$theta[length(df$theta)],df[,test][length(df[,test])],unique(df[,split.par]))
        }
    }
    plot(0,type='n',ylim=range(by.params$egger/by.params$begg,na.rm=TRUE),xlim=range(thetas),xlab=expression(theta),ylab='egger power / begg power',main='ratio')
    for(df in split(by.params,by.params[,split.par])) {
        ratios <- df$egger / df$begg
        lines(df$theta, ratios)
        text(thetas[which.max(ratios)],max(ratios),unique(df[,split.par]))
    }
    if(!is.null(filename)) dev.off()
    par(op)
}

## 5a. S pareto
set.seed(0)
alpha <- .1
## beta <- .8
B <- 1e4
n <- 20
Z.distr <- Z.normal
shapes <- seq(2.01,5,len=3) |> round(digits=2)
thetas <- seq(0,1.5,len=20)
params <- expand.grid(theta=thetas,shape=shapes,stringsAsFactors=FALSE)
res.power <- sapply(1:nrow(params),function(i)sim.power(Z.distr=Z.distr,par.S=S.pareto(c(location=1,shape=params[i,'shape'])),theta=params[i,'theta'],alpha=alpha, B=B,n=n))
by.params.Z.normal <- cbind(params,t(res.power))
sim.plot(by.params.Z.normal,split.par='shape',thetas,filename='S_pareto.png') 

## 5b S beta
set.seed(0)
Z.distr <- Z.normal
shapes <- seq(1,20,len=3) |> round(digits=2)
thetas <- seq(0,5,len=20)
params <- expand.grid(theta=thetas,shape=shapes,stringsAsFactors=FALSE)
res.power <- sapply(1:nrow(params),function(i)sim.power(Z.distr=Z.distr,par.S=S.beta(c(a=1/params[i,'shape'],b=params[i,'shape'],shift=1)),theta=params[i,'theta'],alpha=alpha, B=B,n=n))
by.params.S.beta <- cbind(params,t(res.power))
sim.plot(by.params.S.beta,split.par='shape',thetas,filename='S_beta.png') 

## 5c Z Student's
set.seed(0)
par.S <- S.unif(c(0,1))
thetas <- seq(0,1,len=20)
dfs <- seq(2.01,3,len=3) |> round(digits=2)
params <- expand.grid(theta=thetas,df=dfs,stringsAsFactors=FALSE)
res.power <- sapply(1:nrow(params),function(i)sim.power(Z.distr=Z.t(params[i,'df']),par.S=par.S,theta=params[i,'theta'],alpha=alpha, B=B,n=n)
)
by.params.Z.t <- cbind(params,t(res.power))
sim.plot(by.params.Z.t,split.par='df',thetas,filename='Z_t.png') 

## 5d Z beta
set.seed(0)
par.S <- S.unif(c(0,1))
thetas <- seq(0,2,len=20)
shapes <- c(.51,.65,.75)# |> round(digits=2)
params <- expand.grid(theta=thetas,shape=shapes,stringsAsFactors=FALSE)
res.power <- sapply(1:nrow(params),function(i)sim.power(Z.distr=Z.beta(a=params[i,'shape'],b=1/(params[i,'shape']-.5)),par.S=par.S,theta=params[i,'theta'],alpha=alpha, B=B,n=n)
)
by.params.Z.beta <- cbind(params,t(res.power))
sim.plot(by.params.Z.beta,split.par='shape',thetas,filename='Z_beta.png') 

## 5e S uniform
set.seed(0)
Z.distr <- Z.normal
shifts <- seq(0,5,len=3)
thetas <- seq(0,1.5,len=20)
params <- expand.grid(theta=thetas,shift=shifts,stringsAsFactors=FALSE)
res.power <- sapply(1:nrow(params),function(i)sim.power(Z.distr=Z.distr,par.S=S.unif(c(0,1)+params[i,'shift']),theta=params[i,'theta'],alpha=alpha, B=B,n=n))
by.params.S.uniform <- cbind(params,t(res.power))
sim.plot(by.params.S.uniform,split.par='shift',thetas,filename='S_unif.png') 
