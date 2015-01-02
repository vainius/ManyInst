library(foreach)
# library(Rcpp)
library(ggplot2)
library(glmnet)
library(doParallel)
library(plyr)

GenerateCorrVars <- function(n.obs, n.vars, c){ 
  #gaussian x with cov_ij = c^|j - k|
  
  sigma.m <- foreach(i = 1:n.vars, .combine = rbind)%do%{
    c^abs((seq(1-i, n.vars-i, by = 1)))
  }
  sigma.chol <- chol(sigma.m)
  
  x.uncorr <- matrix(rnorm(n.obs*n.vars), nrow = n.obs, ncol = n.vars)
  x.corr <-  tcrossprod(x.uncorr, sigma.chol)
  x.corr  
}

RJIVEFit <- function(x, Z, L = NULL, par = F){
  if(is.null(L)) L <- sd(x) * sqrt(ncol(Z)) * diag(ncol(Z))
  
  tZZ <- crossprod(Z)
  tLL <- crossprod(L)
  tZX <- crossprod(Z, x)
  
  if(par == F){ 
    OptInst <- rep(NA, length(x))
    for(i in 1:length(x)){
      OptInst[i] <- (t(Z[i, ]) %*% chol2inv(chol(tZZ - tcrossprod(Z[i, ], Z[i, ]) + tLL))) %*% (tZX - tcrossprod(Z[i, ], x[i, ]))     
    }
  } else{
    OptInst <- foreach(i = 1:length(x), .combine = c)%dopar%{
      (t(Z[i, ]) %*% chol2inv(chol(tZZ - tcrossprod(Z[i, ], Z[i, ]) + tLL))) %*% (tZX - tcrossprod(Z[i, ], x[i, ]))     
    }
  }
  
  OptInst
}


RJIVEFitGlmnet <- function(x, Z, nfolds = 20, par = T, lamb.heur = F){
  #in original RJIVE nfolds = n 
  cv.obj <- cv.glmnet(x = Z, y = x, nfolds = nfolds, keep = T, parallel = par, alpha = 0)
  lamb.nr <- max(2, which(cv.obj$lambda == cv.obj$lambda.1se)) #sometimes lambda 1se is too harsh
  if(lamb.heur == T) lamb.nr <- which.min(abs(sd(x) * sqrt(ncol(Z)) - cv.obj$lambda))
  
  lamb.reverse <- ncol(cv.obj$fit.preval) - lamb.nr + 1
  cv.obj$fit.preval[, lamb.reverse]
}


CoefOptInst <- function(y, x, i){
  solve(crossprod(i, x)) %*% crossprod(i, y)
}


OrtoProj <- function(i, W){
  #can put L2 here
  if(class(W) == "numeric"){
    i.orto <- i - W %*% (solve(crossprod(W)) %*% crossprod(W, i))
    i.orto
  } else {
    if(ncol(W) > 0){
      i.orto <- i - W %*% (chol2inv(chol(crossprod(W))) %*% crossprod(W, i))
      i.orto
    } else {
      i
    }
  }
}

# CoefOptInst(y, x, OrtoProj(i, W))
# maybe some corrections needed for many exog

GenerateInst <- function(n.obs, n.inst, W, eps.prop = 0.5, alpha = 1){
  # Z = Wb3 + eps_z 
  signal <- 1/seq(1, ncol(W))^alpha
  
  b3 <- signal * abs(matrix(rnorm(ncol(W) * n.inst), nrow = ncol(W), ncol = n.inst))
  eps.z <- matrix(rnorm(n.obs * n.inst, sd = sqrt(eps.prop)), nrow = n.obs, ncol = n.inst)
  
  b3 <- sqrt(1 - eps.prop) * b3/sqrt(mean(apply(W%*%b3, 2, var)))
  
  Z <- W%*%b3 + eps.z
  Z
}

GenerateX <- function(Z, W, Zg1.prop = 0.25, eps, eps.prop = 0.5, alpha.z = 1, alpha.w = 1){
  # x = Zg1 + Wb2 + eps_x 
  eps.prop <- eps.prop
  
  g1 <- 1/seq(1, ncol(Z))^alpha.z
  b2 <- 1/seq(1, ncol(W))^alpha.w
  
  g1 <- sqrt(Zg1.prop) * g1/sd(Z%*%g1)
  b2 <- sqrt(1 - Zg1.prop - eps.prop) * b2/sd(W%*%b2)
  
  x <- Z%*%g1 + W%*%b2 + eps
  x
}

GenerateY <- function(x, W, W.prop = 0.25, eps, alpha = 1){
  # y = xd + Wb1 + eps_y
  
  d <- 1
  b1 <- 1/seq(1, ncol(W))^alpha
  
  other.var <- var(x*d + eps)
  W.var <- W.prop/(1 - W.prop)*other.var
  
  b1 <- sqrt(W.var) * b1/sd(W%*%b1)
  
  y <- x*d + W%*%b1 + eps
  y
}

GetLassoFit <- function(X, Y, par = F, nfol = 10){
  cv.obj <- cv.glmnet(x = X, y = Y, intercept = F, parallel = par, nfolds = nfol)
  
  lamb.nr <- max(2, which(cv.obj$lambda == cv.obj$lambda.1se)) #sometimes lambda 1se is too harsh
  res <- predict(cv.obj, newx = X, s = cv.obj$lambda[lamb.nr])
  as.numeric(res)
}

GetLassoIDs <- function(X, Y, keepX = NULL){
  
  nonzero.ids <- list()
  
  if(is.null(keepX)){
    
    cv.obj <- cv.glmnet(x = X, y = Y, intercept = F)
    nonzero.ids[["lambda.1se"]] <- which(as.numeric(coef(cv.obj))[-1] != 0) #first coef is intercept
    
    res <- (Y - predict(cv.obj, newx = X))
    lamb.opt <- 2 * sd(res) * sqrt(2 * log(ncol(X) * nrow(X)) / nrow(X))
    lamb.closest <- cv.obj$lambda[which.min(abs(lamb.opt - cv.obj$lambda))]
    nonzero.ids[["rate.optimal"]] <- which(as.numeric(coef(cv.obj, s = lamb.closest))[-1] != 0)
    
  } else {
    Xplus <- cbind(keepX, X)
    cv.obj <- cv.glmnet(x = Xplus, y = Y, intercept = F)
    nonzero.ids[["lambda.1se"]] <- which(as.numeric(coef(cv.obj))[-c(1:(1 + ncol(as.matrix(keepX))))] != 0)
    
    res <- (Y - predict(cv.obj, newx = Xplus))
    lamb.opt <- 2 * sd(res) * sqrt(2 * log(ncol(Xplus) * nrow(Xplus)) / nrow(Xplus))
    lamb.closest <- cv.obj$lambda[which.min(abs(lamb.opt - cv.obj$lambda))]
    nonzero.ids[["rate.optimal"]] <- which(as.numeric(coef(cv.obj, s = lamb.closest))[-1] != 0)
  }
  
  nonzero.ids
}

OrtoRJIVEResults <- function(elim.ids, method.name = "y.lasso.1se", y = y, x = x, W = W, x.fit = x.fit){
  if (length(elim.ids) >= length(y)) {
    #could also choose (n-1) elim.ids at random     
    res.vect <- 0
  } else {
    res.vect <- CoefOptInst(y, x, OrtoProj(i = x.fit, W[, elim.ids]))
  }
  data.frame(coef = res.vect, n_elim = length(elim.ids), method = method.name, elim_ids = paste(elim.ids, collapse = ";"))
}


GetRidgeLambda <- function(X, Y){
  cv.obj <- cv.glmnet(x = X, y = Y, intercept = F, alpha = 0)
  
  lamb.nr <- max(2, which(cv.obj$lambda == cv.obj$lambda.1se)) #sometimes lambda 1se is too harsh
  cv.obj$lambda[lamb.nr]
}


GetRidgeFit <- function(X, Y){
  cv.obj <- cv.glmnet(x = X, y = Y, intercept = F, alpha = 0)
  
  lamb.nr <- max(2, which(cv.obj$lambda == cv.obj$lambda.1se)) #sometimes lambda 1se is too harsh
  res <- predict(cv.obj, newx = X, s = cv.obj$lambda[lamb.nr])
  as.numeric(res)
}


OrtoRJIVEResultsRidge <- function(W.ridge = W.ridge, method.name = "y.ridge.1se", y = y, x = x, x.fit = x.fit){
  res.vect <- CoefOptInst(y, x, OrtoProj(i = x.fit, W.ridge))
  data.frame(coef = res.vect, n_elim = NA, method = method.name, elim_ids = NA)
}


makeOrtoInst <- function(Z, W){
  
  Zort.list <- list()
  
  for(i in 1:ncol(Z)){
    Zi <- Z[, i]
    Zi.fit <- GetLassoFit(X = W, Y = Zi, par = F, nfol = 4)
    Zi.ort <- OrtoProj(i = Zi, W = Zi.fit)
    Zort.list[[i]] <- Zi.ort
  }
  
  do.call(cbind, Zort.list)
}


LassoElim <- function(y = y, x = x, W = W, Z = Z, x.fit = x.fit, keep = F){
  
  #Eliminate W's using LASSO
  # y = xd + Wb1 + eps_y
  # x = Zg1 + Wb2 + eps_x
  # Z = Wb3 + eps_z
  
  res.tmp <- list()
  
  if (keep == F) {
    ids.lasso.elim <- GetLassoIDs(X = W, Y = y)
  } else {
    ids.lasso.elim <- GetLassoIDs(X = W, Y = y, keepX = x)
  }
  
  ids.1se.y <- ids.lasso.elim[["lambda.1se"]]
  ids.rate.y <- ids.lasso.elim[["rate.optimal"]]
  
  if (keep == F) {
    ids.lasso.elim <- GetLassoIDs(X = W, Y = x)
  } else {
    ids.lasso.elim <- GetLassoIDs(X = W, Y = x, keepX = x.fit)
  }
  ids.1se.x <- ids.lasso.elim[["lambda.1se"]]
  ids.rate.x <- ids.lasso.elim[["rate.optimal"]]
  
  if (keep == F) { 
    ids.lasso.elim <- GetLassoIDs(X = W, Y = x.fit)
  } else {
    ids.lasso.elim <- GetLassoIDs(X = W, Y = x.fit, keepX = Z)
  }
  ids.1se.xfit <- ids.lasso.elim[["lambda.1se"]]
  ids.rate.xfit <- ids.lasso.elim[["rate.optimal"]]
  
  ids.1se.yx <- union(ids.1se.y, ids.1se.x)
  ids.rate.yx <- union(ids.rate.y, ids.rate.x)
  ids.1se.yxxfit <- union(union(ids.1se.y, ids.1se.x), ids.1se.xfit)
  ids.rate.yxxfit <- union(union(ids.rate.y, ids.rate.x), ids.rate.xfit)
  
  res.tmp[["y.lasso.1se"]] <- OrtoRJIVEResults(ids.1se.y, "y.lasso.1se", y = y, x = x, W = W, x.fit = x.fit)      
  res.tmp[["y.lasso.rateopt"]] <- OrtoRJIVEResults(ids.rate.y, "y.lasso.rateopt", y = y, x = x, W = W, x.fit = x.fit)
  res.tmp[["x.lasso.1se"]] <- OrtoRJIVEResults(ids.1se.x, "x.lasso.1se", y = y, x = x, W = W, x.fit = x.fit)      
  res.tmp[["x.lasso.rateopt"]] <- OrtoRJIVEResults(ids.rate.x, "x.lasso.rateopt", y = y, x = x, W = W, x.fit = x.fit)
  res.tmp[["xfit.lasso.1se"]] <- OrtoRJIVEResults(ids.1se.xfit, "xfit.lasso.1se", y = y, x = x, W = W, x.fit = x.fit)      
  res.tmp[["xfit.lasso.rateopt"]] <- OrtoRJIVEResults(ids.rate.xfit, "xfit.lasso.rateopt", y = y, x = x, W = W, x.fit = x.fit)
  res.tmp[["yx.lasso.1se"]] <- OrtoRJIVEResults(ids.1se.yx, "yx.lasso.1se", y = y, x = x, W = W, x.fit = x.fit)      
  res.tmp[["yx.lasso.rateopt"]] <- OrtoRJIVEResults(ids.rate.yx, "yx.lasso.rateopt", y = y, x = x, W = W, x.fit = x.fit)
  res.tmp[["yxxfit.lasso.1se"]] <- OrtoRJIVEResults(ids.1se.yxxfit, "yxxfit.lasso.1se", y = y, x = x, W = W, x.fit = x.fit)      
  res.tmp[["yxxfit.lasso.rateopt"]] <- OrtoRJIVEResults(ids.rate.yxxfit, "yxxfit.lasso.rateopt", y = y, x = x, W = W, x.fit = x.fit)
  
  res.bind <- do.call(rbind, res.tmp)
  if(keep == T){
    res.bind$method <- paste0(res.bind$method, ".keep")
  }
  res.bind
}


CalcConcPar <- function(x, opt.inst = x.fit){
  inst.fit <- chol2inv(chol(crossprod(opt.inst))) %*% crossprod(opt.inst, x)
  res <- x - inst.fit
  conc.par <- length(x) * var(inst.fit)/var(res)
  conc.par
}


CalcConcPar <- function(x, opt.inst = x.fit){
  inst.fit <- chol2inv(chol(crossprod(opt.inst))) %*% crossprod(opt.inst, x)
  res <- x - inst.fit
  conc.par <- length(x) * var(inst.fit)/var(res)
  conc.par
}

RmOutliers <- function(vect, mad = F){
  
  MAD <- median(abs(vect - median(vect)))
  Up <- median(vect) + 5 * MAD
  Down <- median(vect) - 5 * MAD
  
  vect[vect > Down & vect < Up]
}

#
#Optimized test-functions
# 

# RJIVEFitFast <- function(x, Z, L = NULL, par = F){
#   if(is.null(L)) L <- sd(x) * sqrt(ncol(Z)) * diag(ncol(Z))
#   
#   tZZ <- crossprod(Z)
#   tLL <- crossprod(L)
#   tZX <- crossprod(Z, x)
#   tZZLL.chol <- chol((tZZ + tLL))
#   
#   if(par == F){ 
#     OptInst <- rep(NA, length(x))
#     for(i in 1:length(x)){
#       OptInst[i] <- (t(Z[i, ]) %*% chol2inv(cholDowndate(tZZLL.chol, Z[i, ]))) %*% (tZX - tcrossprod(Z[i, ], x[i, ]))     
#     }
#   } else{
#     OptInst <- foreach(i = 1:length(x), .export = 'cholDowndate', .combine = c)%dopar%{
#       (t(Z[i, ]) %*% chol2inv(cholDowndate(tZZLL.chol, Z[i, ]))) %*% (tZX - tcrossprod(Z[i, ], x[i, ]))     
#     }
#   }
#   
#   OptInst
# }
# 
# RJIVEFitCpp <- function(x, Z, L = NULL, par = F){
#   if(is.null(L)) L <- sd(x) * sqrt(ncol(Z)) * diag(ncol(Z))
#   
#   tZZ <- crossprod(Z)
#   tLL <- crossprod(L)
#   tZX <- crossprod(Z, x)
#   tZZLL.chol <- as.matrix(chol((tZZ + tLL)))
#   
#   if(par == F){ 
#     OptInst <- rep(NA, length(x))
#     for(i in 1:length(x)){
#       OptInst[i] <- (t(Z[i, ]) %*% chol2inv(cholDowndateCpp(tZZLL.chol, Z[i, ]))) %*% (tZX - tcrossprod(Z[i, ], x[i, ]))     
#     }
#   } else{
#     OptInst <- foreach(i = 1:length(x), .packages = 'Rcpp', .export = 'cholDowndateCpp', .combine = c)%dopar%{
#       (t(Z[i, ]) %*% chol2inv(cholDowndateCpp(tZZLL.chol, Z[i, ]))) %*% (tZX - tcrossprod(Z[i, ], x[i, ]))     
#     }
#   }
#   
#   OptInst
# }
# 
# 
# cholDowndate <- function(R, x){
#   p <- length(x)
#   
#   for (k in 1:(p-1)){
#     r <- sqrt(R[k, k]^2 - x[k]^2)
#     c <- r / R[k, k]
#     s <- x[k] / R[k, k]
#     R[k, k] <- r
#     R[k, c((k+1):p)] <- (R[k, c((k+1):p)] - s * x[(k+1):p]) / c
#     x[(k+1):p] <- c * x[(k+1):p] - s * R[k, (k+1):p]
#   }
#   
#   R[k, k] <- sqrt(R[p, p]^2 - x[p]^2)
#   
#   R
# }
# 
# 
# cppFunction('NumericMatrix cholDowndateCpp(NumericMatrix Rinput, NumericVector xinput) {
#             
#             NumericMatrix R = Rinput;         
#             NumericVector x = xinput;
#             int p = x.size();
#             
#             for (int k = 1; k < p; k++) {
#             double r = sqrt(pow(R(k, k), 2.0) - pow(x[k], 2.0));
#             double con = r / R(k, k);
#             double s = x[k] / R(k, k);
#             R(k, k) = r;
#             
#             for (int z = k+1; z <= p; z++) {
#             R(k, z) = R(k, z) - s * x[z] / con;
#             x[z] = con * x[z] - s * R(k, z);
#             }
#             }
#             
#             R(p, p) = sqrt(pow(R(p, p), 2.0) - pow(x[p], 2.0));
#             
#             return R;
#             }')