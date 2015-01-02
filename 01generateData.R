source("00functions.R")
cl <- makeCluster(4)
registerDoParallel(cl)


SimulateGrand <- function(n.obs.v = c(100, 300, 600), n.sim = 50, glob.alpha = 2, inst.alpha = 1, Zorto = T){
  #alpha govers sparsity. 1/seq(1:k)^(alpha)  
  #Zorto = T ortogonalize instruments before running RJIVE   
  
  # y = xd + Wb1 + eps_y
  # x = Zg1 + Wb2 + eps_x
  # Z = Wb3 + eps_z
  
  foreach(n = n.obs.v, .combine = rbind)%do%{
    
    print(paste(glob.alpha, n))
    W <- GenerateCorrVars(n.obs = n, n.vars = 2*n, c = 0.5)
    
    foreach(sim.nr = 1:n.sim, .combine = rbind)%dopar%{
      
      source("00functions.R")
      
      Z <- GenerateInst(n.obs = n, n.inst = n, W = W, eps.prop = 0.5,
                        alpha = glob.alpha)
      
      cor.epsxy <- 0.5
      eps.x <- rnorm(n, sd = sqrt(0.5))
      eps.y <- cor.epsxy*eps.x + sqrt((1 - cor.epsxy^2))*rnorm(n, sd = sqrt(0.5))
      
      x <- GenerateX(Z, W, Zg1.prop = 0.25, eps = eps.x, eps.prop = 0.5, 
                     alpha.z = inst.alpha, alpha.w = glob.alpha)
      y <- GenerateY(x, W, W.prop = 0.25, eps = eps.y, alpha = glob.alpha)
      
      #Fit RJIVE       
      if (Zorto == T) {
        Z.orto <- makeOrtoInst(Z = Z, W = W)
        x.fit <- RJIVEFitGlmnet(x = x, Z = Z.orto, par = F, lamb.heur = T)
      } else {
        x.fit <- RJIVEFitGlmnet(x = x, Z = Z, par = F, lamb.heur = T)
      }
      
      #Find d
      res.list <- list()
      
      #Use RJIVE inst and eliminate W's      
      res.vect <- CoefOptInst(y, x, i = x.fit)
      res.list[["no.elim"]] <- data.frame(coef = res.vect, n_elim = 0, 
                                          method = "no.elim", elim_ids = NA)
      
      #Eliminate W's using LASSO
      res.list[["lasso"]] <- LassoElim(y = y, x = x, W = W, Z = Z, x.fit = x.fit)
      
      #Eliminate W's using Ridge
      W.ridge <- GetRidgeFit(X = W, Y = y)
      res.list[["y.ridge"]] <- OrtoRJIVEResultsRidge(W.ridge = W.ridge, method.name = "y.ridge.1se", 
                                                     y = y, x = x, x.fit = x.fit)
      
      
      #combine simulation results
      res.dt <- do.call(rbind, res.list)
      res.dt$n.obs <- n
      res.dt$sim.nr <- sim.nr
      res.dt$alpha <- glob.alpha
      res.dt$alpha.inst <- inst.alpha
      res.dt$Zorto <- as.character(Zorto)
      res.dt
      
    }    
  }
}


#Do the simulations
#Specify alpha's, number of obs., number of sims 

set.seed(1000)
sim.res <- list()

#Same alpha for inst. and exog. variables 

for (alpha in c(0.5, 1, 2)){
  print(alpha)
  sim.res[[paste(alpha, "F")]] <- SimulateGrand(n.obs.v = c(500, 1000), n.sim = 45, 
                                                glob.alpha = alpha, inst.alpha = alpha, Zorto = F)
    sim.res[[paste(alpha, "T")]] <- SimulateGrand(n.obs.v = c(500), n.sim = 15, 
                                                  glob.alpha = alpha, inst.alpha = alpha, Zorto = T)
}

#Instruments dense, exog. variables sparse 

print("0.5, 2")
sim.res[[paste("052", "F")]] <- SimulateGrand(n.obs.v = c(500, 1000), n.sim = 45, 
                                              glob.alpha = 3, inst.alpha = 1, Zorto = F)
sim.res[[paste("052", "T")]] <- SimulateGrand(n.obs.v = c(500), n.sim = 15, 
                                              glob.alpha = 3, inst.alpha = 1, Zorto = T)

print("1, 3")
sim.res[[paste("052", "F")]] <- SimulateGrand(n.obs.v = c(500, 1000), n.sim = 60, 
                                              glob.alpha = 3, inst.alpha = 1, Zorto = F)

# 

final.res <- do.call(rbind, sim.res)
saveRDS(final.res, "output/final.res.rds")

