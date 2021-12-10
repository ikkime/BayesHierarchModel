# PC prior setting
U = 5; alpha = 0.05; phi.u = 0.5; phi.alpha = 0.6

pc.prec =list(prior = "pc.prec", param = c(U, alpha))
phi =list(prior = "pc", param = c(phi.u, phi.alpha))

# Bayesian spatial model

mod <- DTH ~ offset(log(CNT)) + 
  
  # Age terms
  f(AGEC, model="rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) +
  
  # Spatial terms
  f(id, model="bym", graph = H, hyper = list(theta = pc.prec)) + 
  
  f(ID_AGE, model="iid", hyper = list(theta = pc.prec))

le.spatial <- inla(mod, family="poisson", data=dat, 
               num.threads = 32, verbose = TRUE, 
               control.fixed = list(mean = 0, prec = 0.01, mean.intercept = 0, prec.intercept = 0.01),
               control.inla = list(strategy = "gaussian"),
               control.predictor = list(compute = TRUE))

# Bayesian spatio-temporal model 1

mod <- DTH ~ 1 + offset(log(CNT)) + 
  
  f(EDU_RANK, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) +
  
  # Spatial terms
  f(id, model="bym", graph = H, hyper = list(theta = pc.prec, phi = phi)) + 
  
  # Age terms
  f(AGEC, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) + 
  
  # Temporal terms
  f(PERIOD, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) +
  
  f(PERIOD1, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec), group = AGEC1) +
  
  f(PERIOD2, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec), group = EDU1_RANK5_1) +  
  
  f(ID_TIME, model = "iid", hyper = list(theta = pc.prec)) + 
  
  f(ID_AGE, model = "iid", hyper = list(theta = pc.prec))

le.st1 <- inla(mod, family="poisson", data=dat, 
               num.threads = 32, verbose = TRUE, 
               control.fixed = list(mean = 0, prec = 0.01, mean.intercept = 0, prec.intercept = 0.01),
               control.inla = list(strategy = "gaussian"),
               control.predictor = list(compute = TRUE))

# Bayesian spatio-temporal model 2

mod <- DTH ~ 1 + offset(log(CNT)) + 

  # Common terms
  AGEC + PERIOD + 
  
  LOG_TOT_POP + LAGE65_PRO + LOG_CTRB + LOG_MED + LN_BUIS +
  
  # Spatial terms and termporal terms
  f(PERIOD1, model = "rw1", scale.model = TRUE, hyper = list(theta = pc.prec), group = AGEC1) +

  f(id, model="bym", graph = H, hyper = list(theta = pc.prec)) + 
  
  f(id2, PERIOD1, model="bym", graph = H, hyper = list(theta = pc.prec)) + 
  
  f(ID_TIME, model = "iid", hyper = list(theta = pc.prec)) +
  
  f(id3, AGEC1, model = "bym", graph = H, hyper = list(theta = pc.prec)) + 
  
  f(ID_AGE, model="iid", hyper = list(theta = pc.prec)) +
  
  f(ID_AGE_TIME, model = "iid", hyper = list(theta = pc.prec))

le.st2 <- inla(mod, family="poisson", data=dat, 
            num.threads = 32, verbose = TRUE, 
            control.fixed = list(mean = 0, prec = 0.01, mean.intercept = 0, prec.intercept = 0.01),
            control.inla = list(strategy = "gaussian"),
            control.predictor = list(compute = TRUE))

# Bayesian spatio-temporal model 3

mod <- DTH ~ offset(log(CNT)) + 
  
  # Common terms
  1 + PERIOD +
  
  # Age terms
  f(AGEC, model="rw1", scale.model = TRUE, hyper = list(theta = pc.prec)) +
  
  f(AGEC1, PERIOD, model = "rw1", hyper = list(theta = pc.prec)) +
  
  f(PERIOD1, model="rw1", scale.model = TRUE, hyper = list(theta = pc.prec), group = AGEC2) +
  
  # Spatial terms
  f(id, model="bym", graph = H, hyper = list(theta = pc.prec)) + 
  
  f(id2, PERIOD, model="bym", graph = H, hyper = list(theta = pc.prec)) +
  
  f(PERIOD2, model="rw1", scale.model = TRUE, hyper = list(theta = pc.prec), group = id3) +
  
  f(ID_AGE, model="iid", hyper = list(theta = pc.prec)) +
  
  f(ID_AGE_TIME, model = "iid", hyper = list(theta = pc.prec))

le.st3 <- inla(mod, family="poisson", data=dat, 
               num.threads = 32, verbose = TRUE, 
               control.fixed = list(mean = 0, prec = 0.01, mean.intercept = 0, prec.intercept = 0.01),
               control.inla = list(strategy = "gaussian"),
               control.predictor = list(compute = TRUE))

# Generation of posterior predictive distribution
mat.post = matrix(,nrow = length(dat[,1]), ncol = nsim)

for (i in c(1:length(dat[,1]))){
  mat.post[i,] <- inla.rmarginal(1000, le.st3$marginals.fitted.values[[i]])
}
