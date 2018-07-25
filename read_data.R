source("functions.R") ## auxiliary functions

## Single culture data
## Tidy data table with all experiments
pyx.ind <- read.csv2("data_single_Pyxidiculla.csv") %>%
    gather(Pyx1A:Pyx3C, key=id, value=counts) %>%
    mutate(batch=as.factor(substr(id,4,4)))
arc.ind <- read.csv2("data_single_Arcella.csv") %>%
    gather(Arcella1A:Arcella3C, key=id, value=counts) %>%
    mutate(batch=as.factor(substr(id,8,8)))
## Matrices for each experiment
pyx1 <- sel.ind(pyx.ind, 1)
pyx2 <- sel.ind(pyx.ind, 2)
pyx3 <- sel.ind(pyx.ind, 3)
arc1 <- sel.ind(arc.ind, 1)
arc2 <- sel.ind(arc.ind, 2)
arc3 <- sel.ind(arc.ind, 3)

## Competition experiments
## Data
comp <- read.csv2("data_competition.csv")
## Matrices for each experiment
experim1 <- sel.exp(comp, batch = 1)
experim2 <- sel.exp(comp, batch = 2)
experim3 <- sel.exp(comp, batch = 3)
experim4 <- sel.exp(comp, batch = 4)
experim5 <- sel.exp(comp, batch = 5)
experim6 <- sel.exp(comp, batch = 6)
experim7 <- sel.exp(comp, batch = 7)
