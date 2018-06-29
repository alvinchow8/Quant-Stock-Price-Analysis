require(roxygen2)

# ===========================================================
# create_GARCH_obj (asset) creates GARCH objects
#' @param asset - individual asset t x 1 vector
#' @param sigma1 - initial sigma1 to initialize GARCH objects
#' @param sigma2 - initial sigma2 to initialize GARCH objects
#' @param T_out - number of time periods for object to predict
#' @export - returns GARCH objects
# ===========================================================
create_GARCH_obj = function(asset,sigma1,sigma2,T_out){
  t = length(asset)
  stan_input = list (T = t,r = asset, sigma1 = sigma1, T_out = T_out)
  GARCH11_fit = stan(file = 'GARCH11_STANCODE.stan', data = stan_input,
                     iter = 500, chains = 10)
  stan_input = list (T = t,r = asset, sigma1 = sigma1, sigma2 = sigma2, T_out = T_out)
  GARCH22_fit = stan(file = 'GARCH22_STANCODE.stan', data = stan_input, iter = 500, chains = 10)
  return(list (G11_obj = GARCH11_fit, G22_obj = GARCH22_fit))
}

# ===========================================================
# price_to_return (fin_data)
#' @param fin_data - financial data t x n matrix
#' @export - returns matrix of returns
# ===========================================================
price_to_return = function (fin_data){
  r_matrix = rep(0,length(fin_data[1,]))
  T_in = length (fin_data[,1])
  for (t in 2:T_in){
    resid_curr = fin_data[t,] - fin_data[t-1,]
    r_matrix = rbind(r_matrix,resid_curr)
  }
  return (r_matrix)
}

# ===========================================================
# plot_PostPrior (GARCH_object)
#' @param GARCH_obj is a GARCH object
#' @param is a matrix of a sector
#' @param bool_11 is true if obj is a GARCH(1,1), false other wise
#' @param sector is the sector which we are looking at (STR)
#' @export : Plots red curve for Posterior, Blue curve for prior
# ===========================================================

plot_PostPrior = function(obj, real_fin_data, bool_11, sector){
  real_fin_data = as.vector(apply(real_fin_data,2,mean))
  if (bool_11 == TRUE){
    plot(density(extract(obj)$r_out),main = cat(sector,
                                                "R_out Density curves - GARCH(1,1)"), col = 'blue')
    lines(density((real_fin_data - mean(real_fin_data))/sd(real_fin_data)), col = 'red')
    plot(density(extract(obj)$alpha0),main = cat(sector,
                                                 "alpha_0 Density curves - GARCH(1,1)"))
    plot(density(extract(obj)$alpha1),main = cat(sector,
                                                 "alpha_1 Density curves - GARCH(1,1)"))
    plot(density(extract(obj)$beta1),main = cat(sector,
                                                "beta_1 Density curves - GARCH(1,1)"))
  }
  else {
      plot(density(extract(obj)$r_out),main = cat(sector,
                                                  "R_out Density curves - GARCH(2,2)"))
      lines(density((real_fin_data - mean(real_fin_data))/sd(real_fin_data)))
      plot(density(extract(obj)$alpha0),main = cat(sector,
                                                   "alpha_0 Density curves - GARCH(2,2)"))
      plot(density(extract(obj)$alpha1),main = cat(sector,
                                                   "alpha_1 Density curves - GARCH(2,2)"))
      plot(density(extract(obj)$alpha2),main = cat(sector,
                                                   "alpha_2 Density curves - GARCH(2,2)"))
      plot(density(extract(obj)$beta1),main = cat(sector,
                                                  "beta_1 Density curves - GARCH(2,2)"))
      plot(density(extract(obj)$beta2),main = cat(sector,
                                                  "beta_2 Density curves - GARCH(2,2)"))
  }
}

# ================================================================
# parameter_extraction (fin_data) : fin_data is a matrix
#' @param fin_data - matrix (t x n) of return data
#' @param sigma1 - initial sigma1 value to initalize GARCH objects
#' @param sigma2 - initial sigma2 value to initialize GARCH objects
#' @param T_out - number of time periods to predict
#' @export
# ===========================================================
parameter_extraction = function(fin_data,sigma1,sigma2,T_out){
  t = length(fin_data[,1])
  data_fin = rep(0,t)
  for (i in 1:t){
    data_fin[i] = mean(as.numeric(fin_data[i,]))
  }
  stan_input = list (T = t,r = data_fin, sigma1 = sigma1, T_out = T_out)
  GARCH11_fit = stan(file = 'GARCH11_STANCODE.stan', data = stan_input,
                     iter = 500, chains = 10)
  GARCH11_param = extract(GARCH11_fit)
  stan_input = list (T = t,r = data_fin, sigma1 = sigma1,sigma2 = sigma2, T_out = T_out)
  GARCH22_fit = stan(file = 'GARCH22_STANCODE.stan', data = stan_input, iter = 500, chains = 10)
  GARCH22_param = extract(GARCH22_fit)

  data_param = list(G11_a0 = mean(GARCH11_param$alpha0),
                      G11_a1 = mean(GARCH11_param$alpha1),
                      G11_b1 = mean(GARCH11_param$beta1),
                      G22_a0 = mean(GARCH22_param$alpha0),
                      G22_a1 = mean(GARCH22_param$alpha1),
                      G22_a2 = mean(GARCH22_param$alpha2),
                      G22_b1 = mean(GARCH22_param$beta1),
                      G22_b2 = mean(GARCH22_param$beta2))
  return (c(data_param,GARCH22_fit,GARCH11_fit))
}

# ===========================================================
# MCMC_diagnostic (GARCH_object, Bool)
#' @param GARCH_object - GARCH object of either GARCH(1,1) or GARCH(2,2)
#' @param bool_11 - True if GARCH(1,1) object, False if GARCH(2,2)
# ===========================================================
MCMC_diagnostic = function(GARCH_object, bool_11){
  GARCH_param = extract(GARCH_object)
  plot (GARCH_param$alpha0,type = 'l',ylab = "Alpha0",main = "Trace Plot")
  plot (GARCH_param$alpha1,type = 'l',ylab = "Alpha1",main = "Trace Plot")
  plot (GARCH_param$beta1,type = 'l',ylab = "Beta1",main = "Trace Plot")
  acf (GARCH_param$alpha0)
  acf (GARCH_param$alpha1)
  acf (GARCH_param$beta1)
  if (bool_11 == FALSE){
    plot (GARCH_param$alpha2,type = 'l',ylab = "Alpha2",main = "Trace Plot")
    plot (GARCH_param$beta2,type = 'l',ylab = "Beta2",main = "Trace Plot")
    acf (GARCH_param$alpha2)
    acf (GARCH_param$beta2)
  }
}

