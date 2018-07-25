source("fitfunctions.R")
## Function for starting values
st.vals3p <- function(){list(r= rnorm(1,2,0.1),
                           k= rnorm(1,200,20),
                           p=runif(1, 0.9, 1),
                           N=ceiling(cbind(50*10*pi*0.125^2/(90*0.25),pyx3*2))+1
                           ) 
}
fit.P3 <-  fit.L(data = pyx3, start.vals= st.vals3p, model="logistic2.jag",
                 n.iter = 5e6, sampArea=10*pi*0.125^2, export="pyx3")
save(fit.P3, file="pyx3.RData")

