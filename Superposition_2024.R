rm(list=ls())
library(MASS)
library(mvtnorm)
library(data.table)
library(ellipse)
library(parallel)
library(matrixcalc)
library(psych)
library(MCMCpack)
library(simex)
library(Matrix)
library(lubridate)
library(truncnorm)
library(matrixStats)
library(numDeriv)
library(tscount)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(gridExtra)
library(forecast)
library(rstudioapi)
library(matrixcalc)




try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)), silent = TRUE)
#setwd("/Users/jiajiekong/Library/CloudStorage/Dropbox/Robert/Poisson_Time_Series")
source("Correlation_B.R")

set.seed(1234)
Y <- rbinom(1000, 1, 0.3)

lambda_t <- function(par, Y, t){
  exp(par[1] + par[2] * t + par[3] * Y[t])
}

Generate_Suposition_Series <- function(par, N_path, lambda_t, Y, n, nsim){
  
  simulate_data_list = list()
  p = par[4]
  phi = par[5]
  
  for (j in 1:nsim){
    
    y = arima.sim(model = list(ar = phi, sd = sqrt(1-phi^2)), n = n*N_path)
    
    y[pnorm(y) > (1-p)] <- 1
    y[pnorm(y) <= (1-p)] <- 0
    
    M = matrix(y, ncol = n, byrow = TRUE)
    X = rep(0, n)
    
    lambda_tt <- lambda_t(par, Y, 1:n)
    lambda_tt_gen <- rpois(n, lambda_tt/p)
    
    for (k in 1:n){
      X[k] = sum(M[1:lambda_tt_gen[k], k])
    }
    
    simulate_data_list[[j]] = X
  }
  
  return(simulate_data_list)
}


Ber_Corr <- function(phi, n){
  
  lst = phi^seq(0,(n-1),by=1)
  asin(lst)/(2*pi)
  
}

# par = c(1, 0.01, 1, 0.5)
# p = 0.5
# N_path = 200
# n = 50
# nsim = 20
# X = Generate_Suposition_Series(par, N_path, lambda_t, Y, n )
#B_corr = Ber_Corr(0.5, 50)

CovX_F = function(lambda_tt, NN, n, B_corr, p){
  
  FF = outer(1:n, 0:NN, FUN = function(X,Y){ppois(Y,lambda_tt[X]/p)})
  E_min = (1-FF) %*% t(1-FF)
  
  CovX = matrix(rep(0, n*n), nrow = n)
  for (i in 1:n){
    for (j in 1:i){
      CovX[i, j] = E_min[i, j] * B_corr[abs(i-j)+1]
      CovX[j, i] = CovX[i, j]
    }
  }
  
  #print(CovX[1:7, 1:7])
  #print(is.positive.definite(CovX))
  return(CovX)
}

obj_function = function(par, X, Y, NN){
  
  p = 1/(1+exp(-par[4]))
  phi = (1/(1+exp(-par[5])) - 1/2) * 2
  print(paste(par[1:3]))
  print(paste(p, phi))
  
  n = length(X)
  B_corr = Ber_Corr(phi, n)
  lambda_tt = lambda_t(par, Y, 1:n)
  
  CovX = CovX_F(lambda_tt, NN, n, B_corr, p)
  
  if (sum(lambda_tt) < 30 || sum(lambda_tt) > 1e7){
    
    Loss = 10^12
    
  }else{
    
    X_hat = rep(0, n); X_hat[1] = lambda_tt[1]
    var_x = rep(0, n); var_x[1] = lambda_tt[1]
    Loss = (X[1] - lambda_tt[1])^2/lambda_tt[1]
    
    print(is.positive.definite(CovX[1:5, 1:5]))
    
    for (k in 2:n){
      
      L = chol(CovX[1:(k-1), 1:(k-1)])
      print(is.positive.definite(CovX[1:(k), 1:(k)]))
      print(k)
      a_tem = forwardsolve(t(L), CovX[k, 1:(k-1)])
      a = backsolve(L, a_tem)
  
      X_hat[k] = (X[1:(k-1)] - lambda_tt[1:(k-1)]) %*% a + lambda_tt[k]
      var_x[k] = (lambda_tt[k]- t(a) %*% CovX[k, 1:(k-1)])
      Loss = Loss + (X[k] - X_hat[k])^2/var_x[k]
  
    }
    
    # X_hat = rep(0, n); X_hat[1] = lambda_tt[1]
    # for (k in 2:n){
    #   
    #   L = chol(CovX[1:(k-1), 1:(k-1)])
    #   a_tem = forwardsolve(t(L), CovX[k, 1:(k-1)])
    #   a = backsolve(L, a_tem)
    #   
    #   X_hat[k] = (X[1:(k-1)] - lambda_tt[1:(k-1)]) %*% a + lambda_tt[k]
    #   #Loss = Loss + (X[k] - X_hat)^2/(lambda_tt[k]- t(a) %*% CovX[k, 1:(k-1)])
    #   
    # }
    # L = chol(CovX)
    # b = forwardsolve(t(L), X - X_hat)
    # Loss = t(b) %*% b
  }
  
  print(Loss)
  return(Loss)
}


main = function(X, Y, NN, par_int, lambda_t, obj_function){
  
  par_int = c(1, 0.01, 1, 0.5, 0.5)
  par_int[4] = log(par_int[4]/(1-par_int[4]))
  par_int[5] = log((par_int[5]/2 + 1/2)/(1-(par_int[5]/2 + 1/2)))
  MLE = optim(par = par_int, fn = obj_function, X = X, Y = Y, NN = NN, lambda_t,
              method = "Nelder-Mead", control = list(factr = 1e7, fnscale=1))
  
  MLE$par[4] = 1/(1+exp(-MLE$par[4]))
  MLE$par[5] = (1/(1+exp(-MLE$par[5])) - 1/2) * 2
  MLE$AIC = -2*MLE$value + 2 * length(par_int)
  MLE$BIC = -2*MLE$value + log(length(X)) * length(par_int)
  
  MLE
}


simulation_learn = function(par_int, N_path, lambda_t, Y, n, simulation_data){
  
  
  NN = 800 #infinity
  N_path = 800 #infinity number of B_t
  nsim = 100
  n = 100
  par = c(1, 0.01, 1, 0.5, 0.5)
  par_int = c(1, 0.01, 1, 0.5, 0.5)
  
  simulation_data = Generate_Suposition_Series(par, N_path, lambda_t, Y, n, nsim)
  
  system.time({
    mc <- getOption("mc.cores", 9)
    simulation_MLE <- mclapply(mc.cores = mc, simulation_data, FUN = main,
                               Y = Y, NN = NN, par_int = par_int, lambda_t, obj_function)

  })
  
  # system.time({
  #   num_cores <- detectCores()
  #   cl <- makeCluster(num_cores)
  #   simulation_MLE <- parLapply(cl, simulation_data, FUN = main, 
  #                              Y = Y, NN = NN, par_int = par_int, lambda_t, obj_function)
  #   stopCluster(cl)
  #   
  # })
  
  par(mfrow = c(1,5))
  ppp = matrix(rep(0,nsim*5), ncol = 5)
  for (i in 1:5){
    for (j in 1:nsim){
      ppp[j, i] = simulation_MLE[[j]]$par[i]
    }
    boxplot(ppp[,i])
    abline(h=par[i], col = "Red", lty = 2)
  }
  
  return(simulation_MLE)
}


main_2 = function(par, N_path, lambda_t, Y, n, NN, nsim, par_int){
  
  simulation_data = Generate_Suposition_Series(par, N_path, lambda_t, Y, n, nsim)
  simulation_MLE = simulation_learn(par_int, N_path, lambda_t, Y, n, simulation_data)
  
  return(list(simulation_data = simulation_data, simulation_MLE = simulation_MLE))
}


system.time({
  
  kk = list()
  p = 0.5
  phi = 0.5
  NN = c(100, 200, 400)  ## number of simulated sample paths
  N_path = c(100, 200, 800)
  nsim = 200   ## number of sample paths for particle filtering approximation
  #Period = c(10, 50) ## length of 1 period
  n = c(50, 100, 300) #c(1e2, 3e2) ## length of the complete timestamp
  par = c(1, 0.01, 1, 0.5, 0.5)
  par_int = c(1, 0.01, 1, 0.5, 0.5)
  
  for (j in 1:3){
    
    kk[[j]] = main_2(par, N_path[j], lambda_t, Y,  n[j], NN[j], nsim, p, par_int, phi)
    
  }
  
})

ppp = array(rep(NA, 3*nsim*3), dim = c(nsim, 3, 3))

for (j in 1:nsim){
  
  try({
    ppp[j,,1] = kk[[1]]$simulation_MLE[[j]]$par
    ppp[j,,2] = kk[[2]]$simulation_MLE[[j]]$par
    ppp[j,,3] = kk[[3]]$simulation_MLE[[j]]$par
  }, silent = TRUE)
  
}


for (j in 1:3){
  
  print(round(apply(ppp[,,j], 2, mean),5))
  print(round(apply(ppp[,,j], 2, sd),5))
  
}

pdf("Poisson_trend_Superposition.pdf", width=6, height=3, paper='special')
par(mfrow = c(1, 3))
ylabb = c('mu', 'beta1', 'beta2')
for (j in 1:3){
  
  boxplot(ppp[,j,1], ppp[,j,2], ppp[,j,3], names = c(50, 100, 300), main = ylabb[j])
  abline(h=par[j], col = "Red", lty = 2)
}
dev.off()

saveRDS(kk, file = "Superposition_July2023")


####most negative problem!
m = 100
grid_u_lambda = exp(seq(-4.6, 2.3025, length.out = m))
p_grid = seq(0.01, .99, length.out = m)
RS = rep(0, m)

for (j in 1:m){
  
  FF = outer(1:100, 0:500, FUN = function(X,Y){ppois(Y,grid_u_lambda[j]/p_grid[X])})
  E_min = diag((1-FF) %*% t(1-FF))
  RS[j] = min( -E_min*p_grid^2/grid_u_lambda[j])
  
}

pdf("plot_Superposition_Most_negative.pdf", width=5, height=5, paper='special')
par(mfrow=c(1,1))
plot(grid_u_lambda, RS, type = 'l', col = 'red',
     xlab=expression(lambda), ylab='Correlation', 
     main = 'Most negative Poisson correlation possible')
dev.off()







