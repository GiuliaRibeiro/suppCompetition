source("fitfunctions.R")

## Function for starting values
st.vals7 <- function(){list(rA= rnorm(1,2,0.1),
                           rP= rnorm(1,2,0.1),
                           kA= rnorm(1,5000,50),
                           kP= rnorm(1,5000,50),    
                           aPA = rnorm(1,1,0.05),
                           aAP = rnorm(1,1,0.05)
                           )
}
## Fitting the model
fit7 <-  fit.LV(data.list = experim7, start.vals= st.vals7, n.iter = 1e8)

save(fit7, file="compfit7.RData")
