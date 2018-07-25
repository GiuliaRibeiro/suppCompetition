## Functions used to fit state-space models with Jags

library(rjags)
library(R2jags)


## Fits a logistic model for each experiment in isolation
fit.L <- function(data, # name of the matrix with counts for both species, returned by sel.ind
                  start.vals= list(r= 2, k= 5000), # initial values
                  N0 = 50,# of cells at Time zero
                  sampArea=10*pi*0.25^2,
                  totArea=0.25*90, #  total sampled area at each replicate and total area of the flask
                  model="logistic.jag", # name of the jags script file with the model
                  par.rec = c("tot", "lamb", "N", "esp", "r", "k", "p"),
                  ...
                   )
{
    n <-  data
    ## Time interval between observations
    dt <- diff(c(0,as.numeric(colnames(data))))
    ## Number of time intervals and replicates
    nInt <- ncol(data) + 1
    nSam <- nrow(data)
    ## List of data to pass to jags
    lista <- list(n=n, dt=dt, nIntervals=nInt, nSamples=nSam,
                  N0=N0, sampArea=sampArea, totArea=totArea)
    if(class(start.vals)=="function")
        init.vals <- start.vals
    else
        init.vals <- function() start.vals
    ## Runs the fit calling jags
    fit1 <- jags.parallel(data = lista,
                          inits = init.vals,
                          parameters.to.save = par.rec,
                          model.file = model,
                          n.chains=4,
                          #export = c("r", "k"),
                          ...
                          )
    return(fit1)
}


## Fits the model from the list of count matrices created by the function sel.exp (see file functions.R)
fit.LV <- function(data.list, # name of the list of matrices with counts for both species, returned by sel.exp
                   start.vals= list(rA= 2, rP= 2, kA= 5000, kP= 5000, aPA = 1, aAP = 1), # initial values
                   A0 = 50, P0=50,# of cells at Time zero
                   sampArea=10*pi*0.25^2,
                   totArea=0.25*90, #  total sampled area at each replicate and total area of the flask
                   model="competition.jag", # name of the jags file with the model
                   par.rec = c("Arc", "Pyx", "lambA", "lambP", "espA", "espP", "rA", "rP", "kA", "kP", "aPA", "aAP", "pA", "pP"),
                   ...
                   )
{
    Arc <-  data.list[[1]]
    Pyx <- data.list[[2]]
    dt <- diff(c(0,as.numeric(colnames(Arc))))
    nInt <- ncol(Arc) + 1
    nSam <- nrow(Arc)
    lista <- list(nA=Arc, nP=Pyx, dt=dt, nIntervals=nInt, nSamples=nSam,
                  A0=A0, P0=P0, sampArea=sampArea, totArea=totArea)
    if(class(start.vals)=="function")
        init.vals <- start.vals
    else
        init.vals <- function() start.vals    
    ## Runs the model
    fit1 <- jags.parallel(data = lista,
                          inits = init.vals,
                          parameters.to.save = par.rec,
                          model.file = model,
                          n.chains=4,
                          #export = c("rA", "rP", "kA", "kP", "aPA", "aAP"),
                          ...
                          )
    return(fit1)
}
