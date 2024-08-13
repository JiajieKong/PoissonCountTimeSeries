Binomial_1 = function(p){ #p = 0.5 in our case
  return(list(Cn = c(0, pbinom(0:1, size = 1, p)), v = p*(1-p)))
}


##H_k##
H_k = function(k, z){
  
  H = rep(0, k)
  if (k == 1){
    R = 1
  } else if (k == 2){
    R = z
  } else {
    H[1] = 1
    H[2] = z
    
    for (j in 3:k){
      H[j] = z * H[j-1] - (j-2) * H[j-2]
    }
    R = H[j]
  }
  R
}

###g_k##
g_k = function(k, Cn){
  
  Phi_inv_Cn = qnorm(Cn)
  H_n = rep(0, length(Phi_inv_Cn))
  
  for (j in 1:length(Phi_inv_Cn)){H_n[j] = H_k(k, Phi_inv_Cn[j])}
  
  sum(exp(-Phi_inv_Cn^2/2) * H_n, na.rm = TRUE) / (factorial(k) * sqrt(2 * pi))
  
}


l_k_cof = function(k, D1, D2){
  
  factorial(k)*g_k(k, D1$Cn)*g_k(k, D2$Cn)/ (sqrt(D1$v * D2$v))
  
}

###correlation plot###
correlation_plot = function(D1, D2, lst){
  
  m = length(lst)
  
  grid_u = lst
  RS = rep(0, m)
  
  for (j in 1:m){
    
    R = 0
    for (k in 1:60){
      
      R = R + l_k_cof(k, D1, D2) * grid_u[j]^k
      
    }
    RS[j] = R

  }
  RS
}


Ber_Corr <- function(phi, n){
  
  lst = phi^seq(0,(n-1),by=1)
  
  correlation_plot(Binomial_1(.5), Binomial_1(.5), lst)
  
}
