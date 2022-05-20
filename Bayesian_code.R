
library("INLA")

U = 5; alpha = 0.05; 

pc.prec =list(prior = "pc.prec", param = c(U, alpha))

# Spatial model
mod0 <- Death ~ 1 + offset(log(POP)) +
  
  f(id, model="bym", graph = adj.rook, hyper = list(theta = pc.prec)) +
  
  f(AGEC, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) +

  f(ID_AGE, model = "iid", hyper = list(theta = pc.prec))

# ST1 model

mod0 <- Death ~ 1 + offset(log(POP)) +
      
    PERIOD +
      
    f(id, model="bym", graph = adj.rook, hyper = list(theta = pc.prec)) +
      
    f(AGEC, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) +
      
    f(ID_AGE_TIME, model = "iid", hyper = list(theta = pc.prec)) 
    
    le0 <- inla(mod0, family="poisson", data=temp, 
                control.fixed = list(mean = 0, prec = 0.01, mean.intercept = 0, prec.intercept = 0.01),
                control.inla = list(strategy = "gaussian"),
                control.compute = list(return.marginals.predictor=TRUE))

# ST2 model 
mod0 <- Death ~ 1 + offset(log(POP)) + 
      
      PERIOD +
      
      # Age/cohort terms
      f(AGEC, model="rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) +
      
      f(AGEC1, PERIOD, model = "rw1", hyper = list(theta = pc.prec)) +
      
      f(PERIOD1, model="rw1", scale.model = TRUE, hyper = list(theta = pc.prec), group = AGEC2) +
      
      # Spatial terms
      f(id, model="bym", graph = H, hyper = list(theta = pc.prec)) + 
      
      f(id2, PERIOD, model="bym", graph = H, hyper = list(theta = pc.prec)) +
      
      f(PERIOD2, model="rw1", scale.model = TRUE, hyper = list(theta = pc.prec), group = id3) +
      
      f(ID_AGE, model="iid", hyper = list(theta = pc.prec)) +
      
      f(ID_AGE_TIME, model = "iid", hyper = list(theta = pc.prec))
    
# Sampling from the posterior predictive distribution    
mat.marginal = matrix(,nrow = length(temp[,1]), ncol = 1000)
    
for (i in c(1:length(temp[,1]))){
      mat.marginal[i,] <- inla.rmarginal(1000, le0$marginals.fitted.values[[i]])
}
