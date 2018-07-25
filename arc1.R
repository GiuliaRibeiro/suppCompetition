source("fitfunctions.R")
## Function for starting values
st.vals1 <- function(){list(r= rnorm(1,2,0.1),
                           k= rnorm(1,200,20),
                           p=runif(1, 0.9, 1),
                           N=ceiling(cbind(50*10*pi*0.25^2/(90*0.25),arc1*2))+1
                           ) 
}
fit.A1 <-  fit.L(data = arc2, start.vals= st.vals1, model="logistic2.jag",
                 n.iter = 5e6, sampArea=10*pi*0.25^2, export="arc1")
save(fit.A1, file="arc1.RData")
save.image()

