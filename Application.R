current_path = dirname(rstudioapi::getSourceEditorContext()$path)
try(setwd(current_path), silent = TRUE)

no_hitting = read.csv("./data/no_hitters_with_covariates.csv", header = TRUE)
Tropical = read.csv("./data/TropicalCyclonesElNino3.csv", header = TRUE)

X1 = Tropical$Storms
Y1 = Tropical$EINo3

X2 = no_hitting$Nohitter
Y2 = no_hitting[,2:3]/1000

##data_visulization##
pdf("App_Tropical_Visulization.pdf", width=8, height=6, paper='special')
par(mfrow = c(2,1), mar = c(4.2, 4, 2, 2))
plot(Tropical$year, Tropical$Storms, type = 'l', 
     ylab = 'Storm', xlab = 'year')
points(Tropical$year, Tropical$Storms, cex = .5)
plot(Tropical$year, Tropical$EINo3, type = 'l', 
     ylab = 'ElNino3', xlab = 'year')
points(Tropical$year, Tropical$EINo3, cex = .5)
# acf(Tropical$Storms)
# pacf(Tropical$Storms)
dev.off()

pdf("App_NoHitting_Visulization.pdf", width=8, height=8, paper='special')
par(mfrow = c(3,1), mar = c(4.2, 4, 2, 2))
plot(no_hitting$Year, no_hitting$Nohitter, type = 'l', 
     ylab = 'no hitting', xlab = 'year')
points(no_hitting$Year, no_hitting$Nohitter)
plot(no_hitting$Year, no_hitting$Number.of.games.played, type = 'l', 
     ylab = '# of games', xlab = 'year')
points(no_hitting$Year, no_hitting$Number.of.games.played)
plot(no_hitting$Year, no_hitting$Mound.height..inch., type = 'l', 
     ylab = 'Mound height inch', xlab = 'year')
points(no_hitting$Year, no_hitting$Mound.height..inch.)
# acf(no_hitting$Nohitter)
# pacf(no_hitting$Nohitter)
dev.off()

##
lambda_t_tropical <- function(par, Y, t){
  exp(par[1] + par[2] * Y[t])
}

lambda_t_no_hitting <- function(par, Y, t){
  exp(par[1] + par[2] * Y[t,1] + par[3] * Y[t,2])
}

##Particle filtering
model_possion_WN_Tropical = function(par_star, org, TT, h){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(par[1:3])
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    #par[4] = (atan(par_star[4])) / pi * 2
    test_pass = TRUE
  }
  
  p = function(t){
    
    probability = exp(par[1] + par[2] * t + par[3] * h[t])
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = 0
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 
model_possion_AR1_Tropical = function(par_star, org, TT, h){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(par[1:3], tan((par[4] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[4] = (atan(par_star[4])) / pi * 2
    test_pass = TRUE
  }
  
  p = function(t){
    
    probability = exp(par[1] + par[2] * t + par[3] * h[t])
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = par[4]
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 
model_possion_AR2_Tropical = function(par_star, org, TT, h){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(par[1:3], tan((par[4:5] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[4:5] = (atan(par_star[4:5])) / pi * 2
    test_pass = TRUE
  }
  
  p = function(t){
    
    probability = exp(par[1] + par[2] * t + par[3] * h[t])
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi1 = par[4]
  phi2 = par[5]
  
  return(list(par = par, p = p, phi1 = phi1, phi2 = phi2, 
              test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 

model_possion_WN_NoHitting = function(par_star, org, TT, h){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(par[1:3])
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    #par[4] = (atan(par_star[4])) / pi * 2
    test_pass = TRUE
  }
  
  p = function(t){
    
    probability = exp(par[1] + par[2] * h[t,1] + par[3] * h[t,2])
    #probability = exp(par[1] )
    
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = 0
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 
model_possion_AR1_NoHitting = function(par_star, org, TT, h){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(par[1:3], tan((par[4] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[4] = (atan(par_star[4])) / pi * 2
    test_pass = TRUE
  }
  
  p = function(t){
    
    probability = exp(par[1] + par[2] * h[t,1] + par[3] * h[t,2])
    #probability = exp(par[1] )
    
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = par[4]
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 
model_possion_AR1_NoHitting_Re = function(par_star, org, TT, h){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(par[1:2], tan((par[3] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[3] = (atan(par_star[3])) / pi * 2
    test_pass = TRUE
  }
  
  p = function(t){
    
    probability = exp(par[1] + par[2] * h[t] )
    #probability = exp(par[1] )
    
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = par[3]
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 

model_possion_AR2_NoHitting = function(par_star, org, TT, h){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(par[1:3], tan((par[4:5] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[4:5] = (atan(par_star[4:5])) / pi * 2
    test_pass = TRUE
  }
  
  p = function(t){
    
    probability = exp(par[1] + par[2] * h[t,1] + par[3] * h[t,2])
    #probability = exp(par[1] )
    
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi1 = par[4]
  phi2 = par[5]
  
  return(list(par = par, p = p, phi1 = phi1, phi2 = phi2, 
              test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 

obj_function_WN = function(par_star, X, model_select, Rp, org, TT, h){ 
  
  model_output = model_select(par_star, org, TT, h)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  Test_Pass = model_output$test_pass
  lag = model_output$lag
  
  N = dim(Rp)[1]
  Z = matrix(rep(0, N*length(X)), nrow = N)
  w = matrix(rep(0, N*length(X)), nrow = N)
  epsilon = matrix(rep(0, N*length(X)), nrow = N)
  log_w = matrix(rep(0, N*length(X)), nrow = N)
  
  for (t in 1:lag){
    P_C_1 = p(t)$PC_1
    phi_1 = phi
    Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2], mean = 0, sd = 1)
    P_x1 = pnorm(P_C_1[X[t] + 2]) - pnorm(P_C_1[X[t] + 1])
    if (t == 1){
      log_w[,t] = log(P_x1) 
    }else{
      log_w[,t] = log_w[,t-1] + log(P_x1)
    }
  }
  
  for (t in (lag+1):length(X)){
    
    P_C = p(t)$PC_1
    sigma = sqrt(1 - phi^2)
    
    ##step 1, find Z_hat##
    Z_hat = Z[,t-lag] * phi
    
    ##step 2, find epsolon##
    lowbound = (P_C[X[t] + 1] - Z_hat) / sigma
    upbound  = (P_C[X[t] + 2] - Z_hat) / sigma
    epsilon[,t] = qtruncnorm(Rp[,t], a = lowbound, b = upbound, mean = 0, sd = 1)
    
    ##step 3, update Z##
    Z[,t] = Z_hat + sigma * epsilon[,t]
    
    ##step 4, update w##
    w_Z = pnorm(upbound) - pnorm(lowbound)
    # w[t] = w[t - 1] * w_Z
    log_w[,t] = log_w[,t - 1] + log(w_Z)
    
  }
  
  log_wt = log_w[,t]
  
  par = round(par, 3)
  
  FF = - log(N) + logSumExp(log_wt)
  
  print(paste(FF, par[1], par[2], par[3], par[4]))
  
  FF
}
obj_function_AR1 = function(par_star, X, model_select, Rp, org, TT, h){ 
  
  model_output = model_select(par_star, org, TT, h)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  Test_Pass = model_output$test_pass
  lag = model_output$lag
  
  N = dim(Rp)[1]
  Z = matrix(rep(0, N*length(X)), nrow = N)
  w = matrix(rep(0, N*length(X)), nrow = N)
  epsilon = matrix(rep(0, N*length(X)), nrow = N)
  log_w = matrix(rep(0, N*length(X)), nrow = N)
  
  for (t in 1:lag){
    P_C_1 = p(t)$PC_1
    phi_1 = phi
    Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2], mean = 0, sd = 1)
    P_x1 = pnorm(P_C_1[X[t] + 2]) - pnorm(P_C_1[X[t] + 1])
    if (t == 1){
      log_w[,t] = log(P_x1) 
    }else{
      log_w[,t] = log_w[,t-1] + log(P_x1)
    }
  }
  
  for (t in (lag+1):length(X)){
    
    P_C = p(t)$PC_1
    sigma = sqrt(1 - phi^2)
    
    ##step 1, find Z_hat##
    Z_hat = Z[,t-lag] * phi
    
    ##step 2, find epsolon##
    lowbound = (P_C[X[t] + 1] - Z_hat) / sigma
    upbound  = (P_C[X[t] + 2] - Z_hat) / sigma
    epsilon[,t] = qtruncnorm(Rp[,t], a = lowbound, b = upbound, mean = 0, sd = 1)
    
    ##step 3, update Z##
    Z[,t] = Z_hat + sigma * epsilon[,t]
    
    ##step 4, update w##
    w_Z = pnorm(upbound) - pnorm(lowbound)
    # w[t] = w[t - 1] * w_Z
    log_w[,t] = log_w[,t - 1] + log(w_Z)
    
  }
  
  log_wt = log_w[,t]
  
  par = round(par, 3)
  
  FF = - log(N) + logSumExp(log_wt)
  
  print(paste(FF, par[1], par[2], par[3], par[4]))
  
  FF
}
obj_function_AR2 = function(par_star, X, model_select, Rp, org, TT, h){ 
  
  model_output = model_select(par_star, org, TT, h)
  p = model_output$p
  phi1 = model_output$phi1
  phi2 = model_output$phi2
  par = model_output$par
  Test_Pass = model_output$test_pass
  lag = model_output$lag
  
  cor_1 = phi1/(1-phi2)
  N = dim(Rp)[1]
  Z = matrix(rep(0, N*length(X)), nrow = N)
  w = matrix(rep(0, N*length(X)), nrow = N)
  epsilon = matrix(rep(0, N*length(X)), nrow = N)
  log_w = matrix(rep(0, N*length(X)), nrow = N)
  
  if (phi2 < -1 || (phi2 > 1 - phi1) || (phi2 > 1 + phi1)){
    FF = -1e16
    return(FF)
  }
  
  t = 1
  P_C_1 = p(t)$PC_1
  Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2], mean = 0, sd = 1)
  P_x1 = pnorm(P_C_1[X[t] + 2]) - pnorm(P_C_1[X[t] + 1])
  log_w[,t] = log(P_x1) 
  
  t = 2
  P_C_1 = p(t)$PC_1
  Z_hat = Z[,t-1] * cor_1
  sigma = sqrt(1-cor_1^2)
  Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2], mean = Z_hat, sd = sigma)
  P_x1 = pnorm(P_C_1[X[t] + 2]) - pnorm(P_C_1[X[t] + 1])
  log_w[,t] = log_w[,t-1] + log(P_x1)
  
  sigma = sqrt(1 - phi1^2 - phi2^2 - phi1*phi2*cor_1)
  
  for (t in 3:length(X)){
    
    P_C = p(t)$PC_1
    
    ##step 1, find Z_hat##
    Z_hat = Z[,t-1] * phi1 + Z[,t-2] * phi2
    
    ##step 2, find epsolon##
    lowbound = (P_C[X[t] + 1] - Z_hat) / sigma
    upbound  = (P_C[X[t] + 2] - Z_hat) / sigma
    epsilon[,t] = qtruncnorm(Rp[,t], a = lowbound, b = upbound, mean = 0, sd = 1)
    
    ##step 3, update Z##
    Z[,t] = Z_hat + sigma * epsilon[,t]
    
    ##step 4, update w##
    w_Z = pnorm(upbound) - pnorm(lowbound)
    # w[t] = w[t - 1] * w_Z
    log_w[,t] = log_w[,t - 1] + log(w_Z)
    
  }
  
  log_wt = log_w[,t]
  
  par = round(par, 3)
  
  FF = - log(N) + logSumExp(log_wt)
  
  print(paste(FF, par[1], par[2], par[3], par[4], par[5]))
  
  FF
}

main_WN = function(X, model_select, para_int_select, N, TT, h){
  
  Rp = matrix(runif(length(X) * N), nrow = N)
  para_int_select_tran = model_select(para_int_select, org = TRUE, TT)$par_tran
  
  system.time({
    MLE_1 = optim(para_int_select_tran, obj_function_WN, 
                  X = X, model_select = model_select, Rp = Rp, org = FALSE, TT = TT, h = h,
                  method = "Nelder-Mead", control = list(factr = 1e7, fnscale=-1))
    
    MLE_1$par_est = model_select(MLE_1$par, org = FALSE, TT)$par
    
    MLE_1$hessian_2 = hessian(func = obj_function_WN, MLE_1$par_est, org = TRUE,
                              X = X, model_select = model_select, Rp = Rp, TT = TT, h = h,
                              method.args=list(eps = 1e-4, d=0.01))
    # MLE_1$se = sqrt(diag(solve(-MLE_1$hessian_2)))
    
    MLE_1$AIC = -2*MLE_1$value + 2 * length(para_int_select_tran)
    MLE_1$BIC = -2*MLE_1$value + log(length(X)) * length(para_int_select_tran)
  })
  
  MLE_1
}
main_AR1 = function(X, model_select, para_int_select, N, TT, h){
  
  Rp = matrix(runif(length(X) * N), nrow = N)
  para_int_select_tran = model_select(para_int_select, org = TRUE, TT)$par_tran
  
  system.time({
    MLE_1 = optim(para_int_select_tran, obj_function_AR1, 
                  X = X, model_select = model_select, Rp = Rp, org = FALSE, TT = TT, h = h,
                  method = "Nelder-Mead", control = list(factr = 1e7, fnscale=-1))
    
    MLE_1$par_est = model_select(MLE_1$par, org = FALSE, TT)$par
    
    MLE_1$hessian_2 = hessian(func = obj_function_AR1, MLE_1$par_est, org = TRUE,
                              X = X, model_select = model_select, Rp = Rp, TT = TT, h = h,
                              method.args=list(eps = 1e-4, d=0.01))
    # MLE_1$se = sqrt(diag(solve(-MLE_1$hessian_2)))
    
    MLE_1$AIC = -2*MLE_1$value + 2 * length(para_int_select_tran)
    MLE_1$BIC = -2*MLE_1$value + log(length(X)) * length(para_int_select_tran)
  })
  
  MLE_1
}
main_AR2 = function(X, model_select, para_int_select, N, TT, h){
  
  Rp = matrix(runif(length(X) * N), nrow = N)
  para_int_select_tran = model_select(para_int_select, org = TRUE, TT)$par_tran
  
  system.time({
    MLE_1 = optim(para_int_select_tran, obj_function_AR2, 
                  X = X, model_select = model_select, Rp = Rp, org = FALSE, TT = TT, h = h,
                  method = "Nelder-Mead", control = list(factr = 1e7, fnscale=-1))
    
    MLE_1$par_est = model_select(MLE_1$par, org = FALSE, TT)$par
    
    MLE_1$hessian_2 = hessian(func = obj_function_AR2, MLE_1$par_est, org = TRUE,
                              X = X, model_select = model_select, Rp = Rp, TT = TT, h = h,
                              method.args=list(eps = 1e-4, d=0.01))
    # MLE_1$se = sqrt(diag(solve(-MLE_1$hessian_2)))
    
    MLE_1$AIC = -2*MLE_1$value + 2 * length(para_int_select_tran)
    MLE_1$BIC = -2*MLE_1$value + log(length(X)) * length(para_int_select_tran)
  })
  
  MLE_1
}

MLE_WN_Tropical = main_WN(X1, model_possion_WN_Tropical, c(2, 0, 0), 500, 1, Y1)
MLE_AR1_Tropical = main_AR1(X1, model_possion_AR1_Tropical, c(2, 0, 0, 0), 500, 1, Y1)
MLE_AR2_Tropical = main_AR2(X1, model_possion_AR2_Tropical, c(2.07, 0.015, -0.284, -0.001, 0), 500, 1, Y1)

MLE_WN_NoHitting = main_WN(X2, model_possion_WN_NoHitting, c(-2, 1, 100), 500, 1, Y2)
MLE_AR1_NoHitting = main_AR1(X2, model_possion_AR1_NoHitting, c(-2, 1, 100, .2), 500, 1, Y2)
MLE_AR2_NoHitting = main_AR2(X2, model_possion_AR2_NoHitting, c(-1.761, 0.896, 80.033, 0.197, 0), 500, 1, Y2)

MLE_AR1_NoHitting_Re = main_AR1(X2, model_possion_AR1_NoHitting_Re, c(-2, 1,.2), 500, 1, Y2$Number.of.games.played)

##residual checking and model checking (PIT)
res_calculation_WN = function(X, par, model_select, org, h){
  
  model_output = model_select(par, org, 1, h)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  
  z_hat_cal = function(X, par, p, phi){
    
    Z_t = rep(0, length(X))
    
    for (t in 1:length(X)){
      
      P_C = p(t)$PC_1
      P_C_2 = p(t)$PC_2
      Z_t[t] = (exp(-P_C[X[t] + 1]^2/2) - exp(-P_C[X[t] + 2]^2/2))/
        (sqrt(2*pi) * (P_C_2[X[t] + 2] - P_C_2[X[t] + 1]))
      
    }
    
    epsilon_t = rep(0, length(X) - 1)
    
    for (t in 2:length(X)){
      
      epsilon_t[t] =  (Z_t[t] ) 
      
    }
    
    return(list(Z_hat = Z_t, epsilon_t_hat = epsilon_t))
  }
  
  z_hat_cal(X, par, p, phi)
  
}
res_calculation_AR1 = function(X, par, model_select, org, h){
  
  
  model_output = model_select(par, org, 1, h)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  
  z_hat_cal = function(X, par, p, phi){
    
    Z_t = rep(0, length(X))
    
    for (t in 1:length(X)){
      
      P_C = p(t)$PC_1
      P_C_2 = p(t)$PC_2
      Z_t[t] = (exp(-P_C[X[t] + 1]^2/2) - exp(-P_C[X[t] + 2]^2/2))/
        (sqrt(2*pi) * (P_C_2[X[t] + 2] - P_C_2[X[t] + 1]))
      
    }
    
    epsilon_t = rep(0, length(X) - 1)
    
    for (t in 2:length(X)){
      
      epsilon_t[t] =  (Z_t[t] - Z_t[t-1] * phi) / sqrt(1 - phi^2)
      
    }
    
    return(list(Z_hat = Z_t, epsilon_t_hat = epsilon_t))
  }
  
  z_hat_cal(X, par, p, phi)
  
}
res_calculation_AR2 = function(X, par, model_select, org, h){
  
  model_output = model_select(par, org, 1, h)
  p = model_output$p
  phi1 = model_output$phi1
  phi2 = model_output$phi2
  par = model_output$par
  cor_1 = phi1/(1-phi2)
  
  z_hat_cal = function(X, par, p, phi){
    
    Z_t = rep(0, length(X))
    
    for (t in 1:length(X)){
      
      P_C = p(t)$PC_1
      P_C_2 = p(t)$PC_2
      Z_t[t] = (exp(-P_C[X[t] + 1]^2/2) - exp(-P_C[X[t] + 2]^2/2))/
        (sqrt(2*pi) * (P_C_2[X[t] + 2] - P_C_2[X[t] + 1]))
      
    }
    
    epsilon_t = rep(0, length(X) - 1)
    
    t = 2
    epsilon_t[2] = Z_t[t] - Z_t[t-1] * cor_1
    sigma = sqrt( 1- phi1^2 - phi2^2 - phi1*phi2*cor_1)
    
    for (t in 3:length(X)){
      
      epsilon_t[t] =  (Z_t[t] - Z_t[t-1] * phi1 - Z_t[t-2]*phi2) / sigma
      
    }
    
    return(list(Z_hat = Z_t, epsilon_t_hat = epsilon_t))
  }
  
  z_hat_cal(X, par, p, phi)
  
}

obj_function_Model_Check_WN = function(par_star, X, model_select, org, h){ #par is vector of 6
  
  model_output = model_select(par_star, org, 1, h)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  
  N = 50
  W_I = rep(0, N)
  Z = rep(0, length(X))
  w = rep(0, length(X))
  epsilon = rep(0, length(X))
  log_w = rep(0, length(X))
  log_wt = rep(0, N)
  
  P_y = array(rep(0, N * length(X) * 100), dim = c(100, length(X), N))
  
  for (i in 1:N){
    
    ##initiate 0##
    P_C_1 = p(1)$PC_1
    phi_1 = phi
    Z[1] = rtruncnorm(1, a = P_C_1[X[1] + 1], b = P_C_1[X[1] + 2], mean = 0, sd = 1)
    w[1] = 1
    P_x1 = pnorm(P_C_1[X[1] + 2]) - pnorm(P_C_1[X[1] + 1])
    log_w[1] = 0
    sigma = 1
    
    Z_hat = 0 ##marginal expectation for the first Obs
    P_y[, 1, i] = pnorm((P_C_1[2:101] - Z_hat)/ 1) - pnorm((P_C_1[1:100] - Z_hat)/ 1)
    
    for (t in 2:length(X)){
      
      P_C = p(t)$PC_1
      
      ##step 1, find Z_hat##
      Z_hat = 0
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma
      epsilon[t] = rtruncnorm(1, a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[t] = Z_hat + sigma * epsilon[t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[t] = log_w[t - 1] + log(w_Z)
      
      ##This part for PIT use##
      for (k in 0:99){
        
        lowbound = (P_C[k + 1] - Z_hat) / sigma
        upbound  = (P_C[k + 2] - Z_hat) / sigma
        P_y[k + 1, t, i] = pnorm(upbound) - pnorm(lowbound)
        
      }
    }
  }
  
  P_y_hat = matrix(rep(0, 100 * length(X)), nrow = 100)
  
  for (j in 1:N){
    
    P_y_hat = P_y_hat + P_y[,,j]
    
  }
  
  P_y_hat = apply(P_y_hat / N, MARGIN = 2, cumsum)
  P_y_hat[100,] = 1
  P_y_hat = rbind(rep(0, length(X)), P_y_hat)
  
  return(P_y_hat = P_y_hat)
}
obj_function_Model_Check_AR1 = function(par_star, X, model_select, org, h){ #par is vector of 6
  
  model_output = model_select(par_star, org, 1, h)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  
  N = 50
  W_I = rep(0, N)
  Z = rep(0, length(X))
  w = rep(0, length(X))
  epsilon = rep(0, length(X))
  log_w = rep(0, length(X))
  log_wt = rep(0, N)
  
  P_y = array(rep(0, N * length(X) * 100), dim = c(100, length(X), N))
  
  for (i in 1:N){
    
    ##initiate 0##
    P_C_1 = p(1)$PC_1
    phi_1 = phi
    Z[1] = rtruncnorm(1, a = P_C_1[X[1] + 1], b = P_C_1[X[1] + 2], mean = 0, sd = 1)
    w[1] = 1
    P_x1 = pnorm(P_C_1[X[1] + 2]) - pnorm(P_C_1[X[1] + 1])
    log_w[1] = 0
    sigma = sqrt(1 - phi^2)
    
    Z_hat = 0 ##marginal expectation for the first Obs
    P_y[, 1, i] = pnorm((P_C_1[2:101] - Z_hat)/ 1) - pnorm((P_C_1[1:100] - Z_hat)/ 1)
    
    for (t in 2:length(X)){
      
      P_C = p(t)$PC_1
      
      ##step 1, find Z_hat##
      Z_hat = Z[t-1] * phi
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma
      epsilon[t] = rtruncnorm(1, a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[t] = Z_hat + sigma * epsilon[t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[t] = log_w[t - 1] + log(w_Z)
      
      ##This part for PIT use##
      for (k in 0:99){
        
        lowbound = (P_C[k + 1] - Z_hat) / sigma
        upbound  = (P_C[k + 2] - Z_hat) / sigma
        P_y[k + 1, t, i] = pnorm(upbound) - pnorm(lowbound)
        
      }
    }
  }
  
  P_y_hat = matrix(rep(0, 100 * length(X)), nrow = 100)
  
  for (j in 1:N){
    
    P_y_hat = P_y_hat + P_y[,,j]
    
  }
  
  P_y_hat = apply(P_y_hat / N, MARGIN = 2, cumsum)
  P_y_hat[100,] = 1
  P_y_hat = rbind(rep(0, length(X)), P_y_hat)
  
  return(P_y_hat = P_y_hat)
}
obj_function_Model_Check_AR2 = function(par_star, X, model_select, org, h){ #par is vector of 6
  
  model_output = model_select(par_star, org, 1, h)
  p = model_output$p
  phi1 = model_output$phi1
  phi2 = model_output$phi2
  par = model_output$par
  
  N = 50
  W_I = rep(0, N)
  Z = rep(0, length(X))
  w = rep(0, length(X))
  epsilon = rep(0, length(X))
  log_w = rep(0, length(X))
  log_wt = rep(0, N)
  
  P_y = array(rep(0, N * length(X) * 100), dim = c(100, length(X), N))
  
  for (i in 1:N){
    
    ##initiate 0##
    P_C_1 = p(1)$PC_1
    phi_1 = phi
    Z[1] = rtruncnorm(1, a = P_C_1[X[1] + 1], b = P_C_1[X[1] + 2], mean = 0, sd = 1)
    w[1] = 1
    P_x1 = pnorm(P_C_1[X[1] + 2]) - pnorm(P_C_1[X[1] + 1])
    log_w[1] = 0
    cor_1 = phi1/(1-phi2)
    sigma = sqrt( 1- phi1^2 - phi2^2 - phi1*phi2*cor_1)
    Z_hat = 0 ##marginal expectation for the first Obs
    P_y[, 1, i] = pnorm((P_C_1[2:101] - Z_hat)/ 1) - pnorm((P_C_1[1:100] - Z_hat)/ 1)
    
    t = 2
    epsilon_t[2] = Z_t[t] - Z_t[t-1] * cor_1
    Z_hat = Z[2] * cor_1
    P_y[, 2, i] = pnorm((P_C_1[2:101] - Z_hat)/ (1-cor_1^2)) - pnorm((P_C_1[1:100] - Z_hat)/ (1-cor_1^2))
    
    for (t in 3:length(X)){
      
      P_C = p(t)$PC_1
      
      ##step 1, find Z_hat##
      Z_hat = Z[t-1] * phi1 + Z[t-2] * phi2
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma
      epsilon[t] = rtruncnorm(1, a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[t] = Z_hat + sigma * epsilon[t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[t] = log_w[t - 1] + log(w_Z)
      
      ##This part for PIT use##
      for (k in 0:99){
        
        lowbound = (P_C[k + 1] - Z_hat) / sigma
        upbound  = (P_C[k + 2] - Z_hat) / sigma
        P_y[k + 1, t, i] = pnorm(upbound) - pnorm(lowbound)
        
      }
    }
  }
  
  P_y_hat = matrix(rep(0, 100 * length(X)), nrow = 100)
  
  for (j in 1:N){
    
    P_y_hat = P_y_hat + P_y[,,j]
    
  }
  
  P_y_hat = apply(P_y_hat / N, MARGIN = 2, cumsum)
  P_y_hat[100,] = 1
  P_y_hat = rbind(rep(0, length(X)), P_y_hat)
  
  return(P_y_hat = P_y_hat)
}

F_bar = function(u, P_y_hat, X){
  
  F_t = function(u, y, t, P_y_hat){
    
    if (u <= P_y_hat[y + 1,t]){
      F_t_u_y = 0
    }else if (u >= P_y_hat[y+2,t]){
      F_t_u_y = 1
    }else{
      F_t_u_y = (u - P_y_hat[y+1, t]) / (P_y_hat[y+2, t] - P_y_hat[y+1, t])
    }
    
    F_t_u_y
  }
  
  kk = 0
  
  for (t in 1:length(X)){
    
    kk = F_t(u, X[t], t, P_y_hat) + kk
    
  }
  
  kk / (length(X))
}#for PIT calculation

Model_Checking = function(MLE, model_select, h, X, f1, f2){
  
  res = f1(X, MLE$par, model_select, org = FALSE, h)
  res$P_y_hat = f2(MLE$par, X, model_select, org = FALSE, h)
  
  return(res)
}

res_WN_Tropical = Model_Checking(MLE_WN_Tropical, model_possion_WN_Tropical, Y1, X1, 
                                 res_calculation_WN, obj_function_Model_Check_WN)
res_AR1_NoHitting = Model_Checking(MLE_AR1_NoHitting_Re, model_possion_AR1_NoHitting_Re, Y2$Number.of.games.played, X2, 
               res_calculation_AR1, obj_function_Model_Check_AR1)

##Report results
plot_res = function(res, X){
  par(mfrow = c(2,2))
  plot(res$epsilon_t_hat, type = 'l', main = 'Residual Series')
  hist(res$epsilon_t_hat, col = 'lightblue', main = 'Residual Histogram')
  #acf(res$epsilon_t_hat)
  pacf(res$epsilon_t_hat)
  grid = seq(0, 1, length.out = 11)
  tt = rep(0, 11)
  for (j in 1:length(grid)){
    tt[j] = F_bar(grid[j], res$P_y_hat, X)
  }
  plot(grid[2:length(grid)], diff(tt), type = "S", ylim = c(0, .2), main = "Model Checking")
  abline(a = .1, b = 0, col = "red", lty = 2)
}

plot_res(res_WN_Tropical, X1)
plot_res(res_AR1_NoHitting, X2)


pdf("App_residual.pdf", width=6, height=5, paper='special')
par(mfrow = c(3,2), mar = c(2, 3, 2, 2))
plot(res_WN_Tropical$epsilon_t_hat, ylab = 'Res Series', main = 'Atlantic Tropical Cyclones',
     type = 'l')
plot(res_AR1_NoHitting$epsilon_t_hat, ylab = 'Res Series', main = 'Baseball No-Hitters',
     type = 'l')
acf(res_WN_Tropical$epsilon_t_hat, ylab = 'Res ACF', main = '')
acf(res_AR1_NoHitting$epsilon_t_hat, ylab = 'Res ACF', main = '')
pacf(res_WN_Tropical$epsilon_t_hat, ylab = 'Res PACF', main = '')
pacf(res_AR1_NoHitting$epsilon_t_hat, ylab = 'Res PACF', main = '')
dev.off()

pdf("App_model_fitness.pdf", width=7, height=3, paper='special')
par(mfrow = c(1,2),mar = c(2, 3, 2, 2) )
grid = seq(0, 1, length.out = 11)
tt = rep(0, 11)
for (j in 1:length(grid)){
  tt[j] = F_bar(grid[j], res_WN_Tropical$P_y_hat, X1)
}
plot(grid[2:length(grid)], diff(tt), type = "S", ylim = c(0, .2), 
     main = 'Atlantic Tropical Cyclones', xlab = '', ylab = '')
abline(a = .1, b = 0, col = "red", lty = 2)

for (j in 1:length(grid)){
  tt[j] = F_bar(grid[j], res_AR1_NoHitting$P_y_hat, X2)
}
plot(grid[2:length(grid)], diff(tt), type = "S", ylim = c(0, .2),  
     main = 'Baseball No-Hitters', xlab = '', ylab = '')
abline(a = .1, b = 0, col = "red", lty = 2)
dev.off()


grid = seq(0, 1, length.out = 11)
tt = rep(0, 11)
for (j in 1:length(grid)){
  tt[j] = F_bar(grid[j], res_WN_Tropical$P_y_hat, X1)
}
model_check_fram = data.frame(gg = grid[2:length(grid)], tt = diff(tt))
p_model_6 = ggplot(model_check_fram) +
  geom_bar(aes(x = gg, y = tt), stat = "identity", fill = "#4983b2") +
  geom_hline(aes(yintercept = .1), linetype = 3) +
  coord_cartesian(ylim = c(0, .20)) +
  labs(x = "PIT-Atlantic Tropical Cyclones", y = "Relative Frequency")

tt2 = rep(0, 11)
for (j in 1:length(grid)){
  tt2[j] = F_bar(grid[j], res_AR1_NoHitting$P_y_hat, X2)
}

model_check_fram = data.frame(gg = grid[2:length(grid)], tt2 = diff(tt))
p_model_2 = ggplot(model_check_fram) +
  geom_bar(aes(x = gg, y = tt), stat = "identity", fill = "#4983b2") +
  geom_hline(aes(yintercept = .1), linetype = 3) +
  coord_cartesian(ylim = c(0, .20)) +
  labs(x = "PIT- Baseball No-Hitters", y = "Relative Frequency")

p_model_2


#uniform test
uniform_test = function(ttt){
  sum(abs(ttt - .1))/10
}

uniform_test_simu = function(tttt){
  #sum(abs(ttt - .1))/10
  cc = rep(0, 10)
  for (j in 1:10){
    cc[j] = sum(((j-1)*0.1 < tttt) *  (tttt< j*0.1))/length(tttt)
  }
  uniform_test(cc)
}

simulate_test = function(L, sim_num, test_stat){
  
  simu_unif = matrix(runif(L * sim_num), nrow = L)
  sum( apply(simu_unif, MARGIN = 2, uniform_test_simu) > test_stat )/sim_num
  
}

simulate_percentile = function(L, sim_num){
  
  simu_unif = matrix(runif(L * sim_num), nrow = L)
  tc = apply(simu_unif, MARGIN = 2, uniform_test_simu) 
  quantile(tc, probs = c(.9, .95, .975, .99, .999))
  
}

uniform_test(diff(tt))
uniform_test(diff(tt2))

simulate_test(53, 1e4, 0.026956)
simulate_test(130, 1e4, 0.02156)

simulate_percentile(53, 1e4)
simulate_percentile(130, 1e4)

pdf("Plot_ModelCheck_Analysis.pdf", width=7, height=3, paper='special')
grid.arrange(arrangeGrob(p_model_6,  p_model_2,ncol = 2), 
             nrow = 1)
dev.off()


c(MLE_WN_NoHitting$AIC, MLE_AR1_NoHitting$AIC, MLE_AR2_NoHitting$AIC)
c(MLE_WN_NoHitting$BIC, MLE_AR1_NoHitting$BIC, MLE_AR2_NoHitting$BIC)

c(MLE_WN_Tropical$AIC, MLE_AR1_Tropical$AIC, MLE_AR2_Tropical$AIC)
c(MLE_WN_Tropical$BIC, MLE_AR1_Tropical$BIC, MLE_AR2_Tropical$BIC)

#sqrt(diag(solve(-MLE_AR2_Tropical$hessian_2)))



