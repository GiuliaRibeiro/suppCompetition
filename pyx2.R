source("fitfunctions.R")
## Function for starting values
st.vals2p <- function(){list(r= rnorm(1,2,0.1),
                           k= rnorm(1,200,20),
                           p=runif(1, 0.9, 1),
                           N=ceiling(cbind(50*10*pi*0.125^2/(90*0.25),pyx2*2))+1
                           ) 
}
fit.P2 <-  fit.L(data = pyx2, start.vals= st.vals2p, model="logistic2.jag",
                 n.iter = 5e6, sampArea=10*pi*0.125^2, export="pyx2")
save(fit.P2, file="pyx2.RData")

