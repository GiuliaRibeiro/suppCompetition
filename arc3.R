source("fitfunctions.R")
## Function for starting values
st.vals3 <- function(){list(r= rnorm(1,2,0.1),
                           k= rnorm(1,200,20),
                           p=runif(1, 0.9, 1),
                           N=ceiling(cbind(50*10*pi*0.25^2/(90*0.25),arc3*2))+1
                           ) 
}
fit.A3 <-  fit.L(data = arc3, start.vals= st.vals3, model="logistic2.jag",
                 n.iter = 5e6, sampArea=10*pi*0.25^2, export="arc3")
save(fit.A3, file="arc3.RData")
save.image()
