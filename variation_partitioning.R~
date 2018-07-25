library(devtools)
#devtools::install_github("mastoffel/rptR", build_vignettes = TRUE)
library(parallel)
library(rptR)
load("cluster/arc1.RData")
load("cluster/arc2.RData")
load("cluster/arc3.RData")
load("cluster/pyx1.RData")
load("cluster/pyx2.RData")
load("cluster/pyx3.RData")
load("cluster/compfit1.RData")
load("cluster/compfit2.RData")
load("cluster/compfit3.RData")
load("cluster/compfit4.RData")
load("cluster/compfit5.RData")
load("cluster/compfit6.RData")
load("cluster/compfit7.RData")
source("functions.R")

## Functions to run parallel simulations
## Function to simulate final population sizes nrep2 times for each of nrep1 draws of posterior parameters
## from fits of single-culture experiments
## i is the index of the experiment (1 to 3)
f1.part <- function(i, sp="arc", nrep1, nrep2, sampArea, totArea=0.25*90, ...){
    dots <- list(...)
    lam <- c()
    if(sp=="arc") sigla <- "A"
    if(sp=="pyx") sigla <- "P"
    for(j in 1:nrep1){
        lam[(j*nrep2-nrep2+1):(j*nrep2)] <- do.call("sim.logist.final",
                                                  args=c(
                                                      list(obj = eval(parse(text=paste("fit.",sigla,i,sep=""))),
                                                           data = eval(parse(text=paste(sp,i, sep=""))),
                                                           nrep = nrep2, sampArea = sampArea, totArea = totArea),
                                                      dots)
                                                  )
    }
    ## convert density to counts
    N <- lam*totArea/sampArea
    data.frame(experiment = i, replicate = rep(1:nrep1, each = nrep2), N )
    }

## Function to simulate final population sizes nrep2 times for each of nrep1 draws of posterior paramenters
## i is the index of the competition experiment (1 to 7)
f2.part <- function(i, nrep1, nrep2, sampArea=10*pi*0.25^2, totArea=0.25*90, ...){
    dots <- list(...)
    lamA <- lamP <- c()
    for(j in 1:nrep1){
        tmp <- do.call("sim.compet.final",
                       args=c(
                           list(obj = eval(parse(text=paste("fit",i,sep=""))),
                                data = eval(parse(text=paste("experim",i, sep=""))),
                                nrep = nrep2, sampArea = sampArea, totArea = totArea),
                           dots)
                       )
        lamA[(j*nrep2-nrep2+1):(j*nrep2)] <- tmp[,1] 
        lamP[(j*nrep2-nrep2+1):(j*nrep2)] <- tmp[,2]
    }
    ## convert densities to counts
    A <- lamA*totArea/sampArea
    P <- lamP*totArea/sampArea
    data.frame(experiment = i, replicate = rep(1:nrep1, each = nrep2), A, P)
    }


## Simulations (parallel)
cl1 <- makePSOCKcluster(4)
clusterExport(cl1, c("logist1", "sim.logist.final",
                     "compet1", "sim.compet.final", 
                     paste("fit",1:7, sep=""),
                     paste("experim",1:7, sep=""),
                     paste("fit.",c("A","P"),rep(1:3,2), sep=""),
                     paste(c("arc","pyx"),rep(1:3,2),sep="")
                     ))

sim.final.arc <- parLapply(cl1, 1:3, fun=f1.part, sp="arc", sampArea=10*pi*0.25^2, nrep1 = 20, nrep2 = 100)
sim.final.pyx <- parLapply(cl1, 1:3, fun=f1.part, sp="pyx",  sampArea=10*pi*0.125^2, nrep1 = 20, nrep2 = 100)
sim.final.comp <- parLapply(cl1, 1:7, fun=f2.part, nrep1 = 20, nrep2 = 100)
stopCluster(cl1)
save.image()

## Analyses ##
## 1. Variance partitioning of final abundances
## Assembling data.frames 
## Arcella
arc.part <- sim.final.arc[[1]]
for(i in 2:length(sim.final.arc))
    arc.part <- rbind(arc.part, sim.final.arc[[i]])
l1 <- nrow(arc.part)
names(arc.part)[3] <- "A"
for(i in 1:length(sim.final.comp))
    arc.part <- rbind(arc.part, sim.final.comp[[i]][,-4])
arc.part$type <- rep(c("single","compet"), c(l1, nrow(arc.part)-l1))
arc.part <- arc.part[,c(4,1:3)]
arc.part$replicate <- with(arc.part, paste(substr(type,1,1),experiment,replicate, sep="."))
arc.part$experiment <- with(arc.part, paste(substr(type,1,1),experiment,sep="."))

## Pyxidiculla
pyx.part <- sim.final.pyx[[1]]
for(i in 2:length(sim.final.pyx))
    pyx.part <- rbind(pyx.part, sim.final.pyx[[i]])
l1 <- nrow(pyx.part)
names(pyx.part)[3] <- "P"
for(i in 1:length(sim.final.comp))
    pyx.part <- rbind(pyx.part, sim.final.comp[[i]][,-3])
pyx.part$type <- rep(c("single","compet"), c(l1, nrow(pyx.part)-l1))
pyx.part <- pyx.part[,c(4,1:3)]
pyx.part$replicate <- with(pyx.part, paste(substr(type,1,1),experiment,replicate, sep="."))
pyx.part$experiment <- with(pyx.part, paste(substr(type,1,1),experiment,sep="."))

## Exploratory boxplots
par(mfrow=c(1,2))
## Expected number of cells in the sample at the end experiments
boxplot(A *(10*pi*0.25^2)/(0.25*90) ~ experiment, data=arc.part, main="Arcella")
boxplot(P * (10*pi*0.25^2)/(0.25*90) ~ experiment, data=pyx.part, main="Pyxidiculla")
par(mfrow=c(1,1))

## Repeatability from Poisson models
## Including mono-specific cultures
## Arcella
rep.A <- rptPoisson(A ~ type + (1 | experiment) + (1 | replicate),
            grname = c("experiment", "replicate", "Fixed", "Residual"),
            data = arc.part, nboot = 200,
            npermut = 0, adjusted = FALSE, parallel=TRUE, ncores = 4)
## Pyxidicolla
rep.P <- rptPoisson(P ~ type + (1 | experiment) + (1 | replicate),
            grname = c("experiment", "replicate", "Fixed", "Residual"),
            data = pyx.part, nboot = 200, 
            npermut = 0, adjusted = FALSE, parallel=TRUE, ncores = 4)
save.image()
summary(rep.P)
summary(rep.A)

## Only for competition experiments
## Arcella
repC.A <- rptPoisson(A ~ 1 + (1 | experiment) + (1 | replicate),
            grname = c("experiment", "replicate", "Residual", "Overdispersion"),
            data = arc.part[arc.part$type=="compet",], nboot = 200,
            npermut = 0, parallel=TRUE, ncores = 4)
## Pyxidicolla
repC.P <- rptPoisson(P ~ 1 + (1 | experiment) + (1 | replicate),
            grname = c("experiment", "replicate", "Residual", "Overdispersion"),
            data = pyx.part[pyx.part$type=="compet",], nboot = 200, 
            npermut = 0, parallel=TRUE, ncores = 4)
save.image()
summary(repC.A)
summary(repC.P)



## 2. Probability of coexistence
## Assembling data.frames 
coex.part <- sim.final.comp[[1]]
for(i in 2:length(sim.final.comp))
    coex.part <- rbind(coex.part, sim.final.comp[[i]])
coex.part$experiment <- with(coex.part, paste("c",experiment,sep="."))
coex.part$replicate <- with(coex.part, paste(experiment,replicate, sep="."))

## Proportions of coexistence
coex.prop <- coex.part %>%
    group_by(experiment, replicate) %>%
    summarise(total=n(),
              n.A = sum(A>0), n.P=sum(P>0),
              n.coex = sum(A>0&P>0), n.nocoex = total - n.coex,
              p.A = n.A/total, p.P= n.P/total, p.coex=n.coex/total) %>%
    as.data.frame()
## Variation in simulated ratios between population sizes
with(coex.part, var(P/A))
## Exploratory plots
boxplot(p.coex~experiment, data=coex.prop)
boxplot(I(n.P/n.A) ~ experiment, data=coex.prop)
boxplot(p.A~experiment, data=coex.prop)
boxplot(p.P~experiment, data=coex.prop)
## Relationship between abundances for each experiment
coex.part %>%
    mutate(repl = factor(substr(replicate,5,5))) %>%
    ggplot(aes(P, A)) +
    geom_point(aes(colour=repl), alpha=0.25) +
    facet_wrap(~experiment)

## Variation partitioning in probability of coexistence
rep.pcoex <- rpt(cbind(n.coex, n.nocoex) ~ 1 + (1 | experiment) + (1 | replicate),
              grname = c("experiment", "replicate", "Residual"),
              datatype = "Proportion",
              data = coex.prop, nboot = 200, 
              npermut = 0, parallel=TRUE, ncores = 4)
rep.pcoexV <- rpt(cbind(n.coex, n.nocoex) ~ 1 + (1 | experiment) + (1 | replicate),
              grname = c("experiment", "replicate", "Residual"),
              datatype = "Proportion", ratio=FALSE,
              data = coex.prop, nboot =0) 
summary(rep.pcoex)
## Variation partitioning in the proportion P/A
rep.pPA <- rpt(cbind(P, A) ~ 1 + (1 | experiment) + (1 | replicate),
              grname = c("experiment", "replicate", "Residual"),
              datatype = "Proportion",
              data = coex.part, nboot = 200, 
              npermut = 0, parallel=TRUE, ncores = 4)
## Variances
rep.pPAV <- rpt(cbind(P, A) ~ 1 + (1 | experiment) + (1 | replicate),
              grname = c("experiment", "replicate", "Residual"),
              datatype = "Proportion", ratio=FALSE,
              data = coex.part, nboot = 200, 
              npermut = 0, parallel=TRUE, ncores = 4)
save.image()
summary(rep.pPA)
summary(rep.pPAV)

## Variance of logits of ratios between observed cell counts 
list(experim1, experim2, experim3, experim4, experim5, experim6, experim7) %>%
    sapply(function(x) c(sum(x[[1]][,ncol(x[[1]])]), sum(x[[2]][,ncol(x[[2]])]))) %>%
    t() %>%
    as.data.frame() %>%
    mutate(total=V1+V2, ratio=V2/total , lratio=log(ratio/(1-ratio))) %>%
    summarize(vratio=var(lratio))
