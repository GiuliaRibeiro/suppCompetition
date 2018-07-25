library(zoo)
library(xts)
library(caTools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rjags)
library(R2jags)

## Creates matrices of counts x time
sel.ind <- function(data, batch.number = 1, time.limit = 20){
    filter(data, Time < time.limit, batch==batch.number) %>%
        mutate(replicate=rep(letters[1:3],length(unique(Time))))%>%
        select(Time, replicate, counts) %>%
        spread(Time, value=counts) %>%
        select(-replicate) %>%
        as.matrix()
}


## Creates a list of matrices of counts x time for each experiment
sel.exp <- function(data = comp, batch.number = 1, time.limit = 20){
    Arc <- filter(data, Time < time.limit, batch==batch.number) %>%
        mutate(replicate=rep(letters[1:3],length(unique(Time))))%>%
        select(Time, replicate, arc) %>%
        spread(Time, value=arc) %>%
        select(-replicate) %>%
        as.matrix()
    Pyx <- filter(data, Time < time.limit&batch==batch.number) %>%
        mutate(replicate=rep(letters[1:3],length(unique(Time))))%>%
        select(Time, replicate, pyx) %>%
        spread(Time, value=pyx) %>%
        select(-replicate) %>%
        as.matrix()
    list(Arc, Pyx)
}

## Calculates the overlap of two posterior distributions that are in a object of class "rjags"
post.overl <- function(obj1, obj2, par.name){
    d1 <- density(obj1$BUGSoutput$sims.list[[par.name]])
    d2 <- density(obj2$BUGSoutput$sims.list[[par.name]])
    f1 <- approxfun(d1, yleft=0, yright=0)
    f2 <- approxfun(d2, yleft=0, yright=0)
    f3 <- Vectorize(function(x) min(f1(x),f2(x)))
    if(min(d1$x)>max(d2$x)|min(d2$x)>max(d1$x))
        x <- 0
    else
        x <- integrate(f3, min(c(d1$x, d2$x)), max(d1$x,d2$x))$value
    return(x)
    }

## Plots posteriors of a given parameter for each experiment
post.plot <- function(lista, par.name, transf=FALSE, legend=TRUE, ...){
    if(transf)
        lista <- lapply(lista, function(x) lapply(x,log))
    lista.dens <- lapply(lista, function(x)density(x[[par.name]]))
    x.minimos <- sapply(lista.dens, function(x)min(x$x))
    x.maximos <- sapply(lista.dens, function(x)max(x$x))
    x.limites <- c(min(x.minimos), max(x.maximos))
    y.minimos <- sapply(lista.dens, function(x)min(x$y))
    y.maximos <- sapply(lista.dens, function(x)max(x$y))
    y.limites <- c(min(y.minimos), max(y.maximos))
    plot(lista.dens[[1]], xlim = x.limites, ylim = y.limites, main=par.name, ...)
    if(length(lista.dens)>1){
        for(i in 2:length(lista.dens))
            lines(lista.dens[[i]], col=i)
        if(legend)
            legend("topright", names(lista.dens), lty=1,
                   col=1:length(lista), bty="n")
    }
}

## Gets the mean and 95% CI of the posterior distributions of the parameters that do not vary in time
dem.pars <- function(object){
    tmp <- object$BUGSoutput$summary
    df <- data.frame(rA=tmp[grep("rA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     rP=tmp[grep("rP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     kA=tmp[grep("kA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     kP=tmp[grep("kP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     aPA=tmp[grep("aPA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     aAP=tmp[grep("aAP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     pA=tmp[grep("^pA",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")],
                     pP=tmp[grep("^pP",rownames(tmp)),][c("mean","sd","2.5%","97.5%","Rhat","n.eff")])
    t(df)
    }
    
## simulates a logistic dynamic as modelled 
logist1 <- function(r, k, p, dt, sampArea, totArea=0.25*90, N0=50, obs=TRUE){
    tot <- N0
    n <- lamb <- c()
    for(t in 1:(length(dt))){
        mu <- tot + (r * tot * (1 - tot/k))*dt[t]
        tot <- rpois(1, mu)
        lamb[t] <- tot*sampArea/totArea
        N <- rpois(1,lamb)
        n[t] <- rbinom(1, N, p)
    }
if(obs)    
    return(n)
else
    return(lamb)
    
}

## sample from posterior of simulated logistic from estimated lambda
post.logist <- function(obj, nrep=1000){
    results <- matrix(ncol=ncol(obj$BUGSoutput$sims.list$lamb), nrow=nrep)
    for(j in 1:nrep){
        i <- sample(1:length(obj$BUGSoutput$sims.list$r), 1)
        lam <- obj$BUGSoutput$sims.list$lamb[i,]
        p <- obj$BUGSoutput$sims.list$p[i]
        N <- rpois(length(lam), lam)
        results[j,] <- rbinom(length(lam), N, p)
    }
    results
}

## Simulates a logistic from their parameters
sim.logist <- function(obj, data, sampArea, totArea=0.25*90, N0=50, obs=FALSE, nrep=1000){
    dt <- diff(c(0,as.numeric(colnames(data))))
    results <- matrix(ncol=length(dt), nrow=nrep)
    for(j in 1:nrep){
        i <- sample(1:length(obj$BUGSoutput$sims.list$r), 1)
        r <- obj$BUGSoutput$sims.list$r[i]
        k <- obj$BUGSoutput$sims.list$k[i]
        p <- obj$BUGSoutput$sims.list$p[i]
        results[j,] <- logist1(r, k, p, dt, sampArea, totArea, N0, obs=obs)
    }
    results
}


## Run simulations of logistic single-species dynamics nrep times with the same parameter values taken from posterior
## Returns only the final population size
sim.logist.final <- function(obj, data, dt, sampArea, totArea=0.25*90, N0=50, obs=FALSE, nrep=100){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data))))
    results <- c()
    i <- sample(1:length(obj$BUGSoutput$sims.list$r), 1)
    for(j in 1:nrep){
        r <- obj$BUGSoutput$sims.list$r[i]
        k <- obj$BUGSoutput$sims.list$k[i]
        p <- obj$BUGSoutput$sims.list$p[i]
        results[j] <- logist1(r, k, p, dt, sampArea, totArea, N0, obs=obs)[length(dt)]
    }
    results
}

## Simulates the competition dynamics with stochasticity
compet1 <- function(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea=0.25*90, A0=50, P0=50, obs=TRUE){
    A <- A0
    P <- P0
    lambA <- lambP <- yA <- yP <- c()
    for(t in 1:(length(dt))){
        muA <- A + (rA * A * (1 - ((A+aPA*P)/kA)))*dt[t]
        muP <- P + (rP * P * (1 - ((P+aAP*A)/kP)))*dt[t]
        A <- ifelse(muA>0, rpois(1, muA),0)
        P <- ifelse(muP>0, rpois(1, muP), 0)
        lambA[t] <- A*sampArea/totArea
        lambP[t] <- P*sampArea/totArea
        nA <- rpois(1,lambA[t])
        nP <- rpois(1,lambP[t])
        yA[t] <- rbinom(1, nA, pA)
        yP[t] <- rbinom(1, nP, pP)
    }
    if(obs)
        cbind(yA,yP)
    else
        cbind(lambA, lambP) 
}

## Simulates the competition dynamics without stochasticity
compet2 <- function(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea=0.25*90, A0=50, P0=50, obs=TRUE){
    A <- P <- c()
    A[1] <- A0
    P[1] <- P0
    for(t in 2:(length(dt))){
        A[t-1] <- ifelse(A[t-1] > 0, A[t-1] , 0)
        P[t-1] <- ifelse(P[t-1] > 0, P[t-1] , 0)
        A[t] <- A[t-1] + (rA * A[t-1] * (1 - ((A[t-1]+aPA*P[t-1])/kA)))*dt[t]
        P[t] <- P[t-1] + (rP * P[t-1] * (1 - ((P[t-1]+aAP*A[t-1])/kP)))*dt[t]        
    }
    lambA <- A*sampArea/totArea
    lambP <- P*sampArea/totArea
    if(obs){
        nA <- rpois(length(A) ,lambA)
        nP <- rpois(length(A) ,lambA)
        yA <- rbinom(length(nA), nA, pA)
        yP <- rbinom(length(nA), nP, pP)
        return(cbind(yA, yP))
    }
    else
        return(cbind(lambA, lambP))
}

## Posterior of competition, from expected posterior of expected population sizes
## Run simulations of competition dynamics, with parameters taken form the posterior distribution
post.compet <- function(obj, nrep=1000){
    nobs <- ncol(obj$BUGSoutput$sims.list$lambA)
    results <- array( dim=c(nobs, 2, nrep))
    for(j in 1:nrep){
        i <- sample(1:nrow(obj$BUGSoutput$sims.list$lambA), 1)
        lambA <- obj$BUGSoutput$sims.list$lambA[i,]
        lambP <- obj$BUGSoutput$sims.list$lambP[i,]
        pA <- obj$BUGSoutput$sims.list$pA[i]
        pP <- obj$BUGSoutput$sims.list$pP[i]
        nA <- rpois(nobs,lambA)
        nP <- rpois(nobs,lambP)
        results[,1,j] <- rbinom(nobs, nA, pA)
        results[,2,j] <- rbinom(nobs, nP, pP)
    }
    results
}

## Run simulations of competition dynamics, with parameters taken from the posterior distribution
sim.compet <- function(obj, data, dt, sampArea=10*pi*0.25^2, totArea=0.25*90, A0=50, P0=50, nrep=1000, obs=FALSE){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data[[1]]))))
    res.a <- matrix(NA, nrep, 3)
    res.est <- array( dim=c(length(dt), 2, nrep))
    res.no.est <- array( dim=c(length(dt), 2, nrep))
    for(j in 1:nrep){
        i <- sample(1:length(obj$BUGSoutput$sims.list$rA), 1)
        rA <- obj$BUGSoutput$sims.list$rA[i]
        kA <- obj$BUGSoutput$sims.list$kA[i]
        pA <- obj$BUGSoutput$sims.list$pA[i]
        rP <- obj$BUGSoutput$sims.list$rP[i]
        kP <- obj$BUGSoutput$sims.list$kP[i]
        pP <- obj$BUGSoutput$sims.list$pP[i]
        aPA <- obj$BUGSoutput$sims.list$aPA[i]
        aAP <- obj$BUGSoutput$sims.list$aAP[i]
        A <- (kA - aAP*kP)/(1-aAP*aPA)
        P <- (kP - aPA*kA)/(1-aAP*aPA)
        estavel <- 1/aPA < kA/kP & kA/kP < aAP
        ## Analytical solution at equilibrium
        res.a[j,] <- c( A, P, estavel)
        ## Simulation with stochasticity at the dt interval
        res.est[,,j] <- compet1(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)
        ## simulacao without stochasticity
        res.no.est[,,j] <- compet2(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)        
    }
    list(analitico=res.a, estocast=res.est, sem.estocast = res.no.est)
}


## Same simulation, but using parameters r and k taken from the pure cultures experiments
sim.compet2 <- function(obj1, obj2, obj3, data, dt, sampArea=10*pi*0.25^2, totArea=0.25*90, A0=50, P0=50, nrep=1000, obs=FALSE){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data[[1]]))))
    res.a <- matrix(NA, nrep, 3)
    res.est <- array( dim=c(length(dt), 2, nrep))
    res.no.est <- array( dim=c(length(dt), 2, nrep))
    for(k in 1:nrep){
        i <- sample(1:length(obj1$BUGSoutput$sims.list$rA), 1)
        j <- sample(1:length(obj2$BUGSoutput$sims.list$r), 1)
        z <- sample(1:length(obj3$BUGSoutput$sims.list$r), 1)
        rA <- obj2$BUGSoutput$sims.list$r[j]
        kA <- obj2$BUGSoutput$sims.list$k[j]
        pA <- obj1$BUGSoutput$sims.list$pA[i]
        rP <- obj3$BUGSoutput$sims.list$r[z]
        kP <- obj3$BUGSoutput$sims.list$k[z]
        pP <- obj1$BUGSoutput$sims.list$pP[i]
        aPA <- obj1$BUGSoutput$sims.list$aPA[i]
        aAP <- obj1$BUGSoutput$sims.list$aAP[i]
        A <- (kA - aAP*kP)/(1-aAP*aPA)
        P <- (kP - aPA*kA)/(1-aAP*aPA)
        estavel <- 1/aPA < kA/kP & kA/kP < aAP
        res.a[k,] <- c( A, P, estavel)
        res.est[,,k] <- compet1(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)
        ## simulacao sem estocasticidade ambiental (para estocasticidade observacional use obs=TRUE)
        res.no.est[,,k] <- compet2(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs) 
    }
    list(analitico=res.a, estocast=res.est, sem.estocast = res.no.est)
}



## Same simulation, but using parameters r and k taken from the pure cultures experiments only for Arcella
sim.compet3 <- function(obj1, obj2, data, dt, sampArea=10*pi*0.25^2, totArea=0.25*90, A0=50, P0=50, nrep=1000, obs=FALSE){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data[[1]]))))
    res.a <- matrix(NA, nrep, 3)
    res.est <- array( dim=c(length(dt), 2, nrep))
    res.no.est <- array( dim=c(length(dt), 2, nrep))
    for(k in 1:nrep){
        i <- sample(1:length(obj1$BUGSoutput$sims.list$rA), 1)
        j <- sample(1:length(obj2$BUGSoutput$sims.list$r), 1)
        rA <- obj2$BUGSoutput$sims.list$r[j]
        kA <- obj2$BUGSoutput$sims.list$k[j]
        pA <- obj1$BUGSoutput$sims.list$pA[i]
        rP <- obj1$BUGSoutput$sims.list$rP[i]
        kP <- obj1$BUGSoutput$sims.list$kP[i]
        pP <- obj1$BUGSoutput$sims.list$pP[i]
        aPA <- obj1$BUGSoutput$sims.list$aPA[i]
        aAP <- obj1$BUGSoutput$sims.list$aAP[i]
        A <- (kA - aAP*kP)/(1-aAP*aPA)
        P <- (kP - aPA*kA)/(1-aAP*aPA)
        estavel <- 1/aPA < kA/kP & kA/kP < aAP
        res.a[k,] <- c( A, P, estavel)
        res.est[,,k] <- compet1(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)
        res.no.est[,,k] <- compet2(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)
    }
    list(analitico=res.a, estocast=res.est, sem.estocast = res.no.est)
}

## Same simulation, but using parameters r and k taken from the pure cultures experiments only for Pixydiculla
sim.compet4 <- function(obj1, obj2, data, dt, sampArea=10*pi*0.25^2, totArea=0.25*90, A0=50, P0=50, nrep=1000, obs=FALSE){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data[[1]]))))
    res.a <- matrix(NA, nrep, 3)
    res.est <- array( dim=c(length(dt), 2, nrep))
    res.no.est <- array( dim=c(length(dt), 2, nrep))
    for(k in 1:nrep){
        i <- sample(1:length(obj1$BUGSoutput$sims.list$rA), 1)
        j <- sample(1:length(obj2$BUGSoutput$sims.list$r), 1)
        rA <- obj1$BUGSoutput$sims.list$rA[i]
        kA <- obj1$BUGSoutput$sims.list$kA[i]
        pA <- obj1$BUGSoutput$sims.list$pA[i]
        rP <- obj2$BUGSoutput$sims.list$r[j]
        kP <- obj2$BUGSoutput$sims.list$k[j]
        pP <- obj1$BUGSoutput$sims.list$pP[i]
        aPA <- obj1$BUGSoutput$sims.list$aPA[i]
        aAP <- obj1$BUGSoutput$sims.list$aAP[i]
        aPA <- obj1$BUGSoutput$sims.list$aPA[i]
        aAP <- obj1$BUGSoutput$sims.list$aAP[i]
        A <- (kA - aAP*kP)/(1-aAP*aPA)
        P <- (kP - aPA*kA)/(1-aAP*aPA)
        estavel <- 1/aPA < kA/kP & kA/kP < aAP
        res.a[k,] <- c( A, P, estavel)
        res.est[,,k] <- compet1(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)
        res.no.est[,,k] <- compet2(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)
    }
    list(analitico=res.a, estocast=res.est, sem.estocast = res.no.est)
}


## Run simulations of competition dynamics nrep times with the same parameter values taken from posterior
## Returns only the final population size
sim.compet.final <- function(obj, data, dt, sampArea=10*pi*0.25^2, totArea=0.25*90, A0=50, P0=50, nrep=100, obs=FALSE){
    if(missing(dt))
        dt <- diff(c(0,as.numeric(colnames(data[[1]]))))
    results <- matrix(NA, nrep, 2)
    i <- sample(1:length(obj$BUGSoutput$sims.list$rA), 1)
    for(j in 1:nrep){
        rA <- obj$BUGSoutput$sims.list$rA[i]
        kA <- obj$BUGSoutput$sims.list$kA[i]
        pA <- obj$BUGSoutput$sims.list$pA[i]
        rP <- obj$BUGSoutput$sims.list$rP[i]
        kP <- obj$BUGSoutput$sims.list$kP[i]
        pP <- obj$BUGSoutput$sims.list$pP[i]
        aPA <- obj$BUGSoutput$sims.list$aPA[i]
        aAP <- obj$BUGSoutput$sims.list$aAP[i]
        results[j,] <- compet1(rA, kA, pA, rP, kP, pP, aPA, aAP, dt, sampArea, totArea, A0, P0, obs=obs)[length(dt),]        
    }
    results
}


## zoo object with predicted and observed values
## For logistic
summary.df <- function(fit.obj, data){
    tmp <- fit.obj$BUGSoutput$summary
    tmp2 <- post.logist(fit.obj)
    zoo(
        data.frame(
            obs1=data[1,],
            obs2=data[2,],
            obs3=data[3,],
            obsm = apply(data,2,mean),
        esp=tmp[grep("lamb",rownames(tmp)),"mean"][-1],
        esp.low=tmp[grep("lamb",rownames(tmp)),"2.5%"][-1],
        esp.up=tmp[grep("lamb",rownames(tmp)),"97.5%"][-1],
        post.low=apply(tmp2,2,quantile,0.025)[-1],
        post.up=apply(tmp2,2,quantile,0.975)[-1]),
        order.by=as.numeric(colnames(data)))
}

## Zoo object with predicted and observed values for competition
summary.df2 <- function(fit.obj, experim, k=1){
    tmp <- fit.obj$BUGSoutput$summary
    tmp2 <- post.compet(fit.obj)
    zoo(
        data.frame(obsA=runmean(apply(experim[[1]],2,mean),k),
                   obsP=runmean(apply(experim[[2]],2,mean),k),
                   espA=tmp[grep("lambA",rownames(tmp)),"mean"][-1],
                   espA.low=tmp[grep("lambA",rownames(tmp)),"2.5%"][-1],
                   espA.up=tmp[grep("lambA",rownames(tmp)),"97.5%"][-1],
                   espP=tmp[grep("lambP",rownames(tmp)),"mean"][-1],
                   espP.low=tmp[grep("lambP",rownames(tmp)),"2.5%"][-1],
                   espP.up=tmp[grep("lambP",rownames(tmp)),"97.5%"][-1],
                   projA.low=apply(tmp2,c(1,2),quantile,0.025, type=9)[,1][-1],
                   projA.up=apply(tmp2,c(1,2),quantile,0.975, type=9)[,1][-1],
                   projP.low=apply(tmp2,c(1,2),quantile,0.025, type=9)[,2][-1],
                   projP.up=apply(tmp2,c(1,2),quantile,0.975, type=9)[,2][-1]),
        order.by=as.numeric(colnames(experim[[1]])))
}

## Only projected trajectories for competition
comp.proj <- function(fit.obj, experim, dt, ...){
    tmp <- sim.compet(fit.obj, dt = dt, data = experim, ...)
    zoo(
        data.frame(
            espA=apply(tmp,c(1,2),mean)[,1][-1],
            espP=apply(tmp,c(1,2),mean)[,2][-1],
            espA.low=apply(tmp,c(1,2),quantile,0.025, type=9)[,1][-1],
            espA.up=apply(tmp,c(1,2),quantile,0.975, type=9)[,1][-1],
            espP.low=apply(tmp,c(1,2),quantile,0.025, type=9)[,2][-1],
            espP.up=apply(tmp,c(1,2),quantile,0.975, type=9)[,2][-1]),
        order.by=cumsum(c(0,dt)))
}

## Mean of lognormal distribution
mu.ln <- function(x,lsd) log(x)-(lsd^2)/2

