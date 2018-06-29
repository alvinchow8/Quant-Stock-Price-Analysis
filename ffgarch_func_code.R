require(rstan)
require(roxygen2)
require(LaplacesDemon)
require(testthat)
options(mc.cores = parallel::detectCores())

#------------------------------------------------------------------------------------------------------------------

#' Helper functionto create returns from prices as input for FFGARCH.Fit()  
#' 
#' @param fin_data Matrix of stock prices
#' @param type Type of financial returns to be calculated (1 = absolute, 2 = relative, 3 = log)
#' @param frame Should function produce returns as data frame (not for model) of L. of Lists

returns_prices <- function(fin_data, type, frame) {
  expect_that(fin_data,is_a("data.frame"))
  expect_that(type,is_a("numeric"))
  expect_that(frame,is_a("logical"))
  expect_equal(type > 0, TRUE)
  expect_equal(type < 4, TRUE)
  
  tot <- nrow(fin_data) - 1
  stocks <- ncol(fin_data) 
  ret_frame <- fin_data
  if (type == 1) {
    for (t in 1:tot) {
      ret_frame[t,] = fin_data[t+1,] - fin_data[t,]
    }
    ret_frame <- ret_frame[1:tot,]
  }
  else if (Type == 2) {
    for (t in 1:tot) {
      ret_frame[t,] = (fin_data[t+1,] - fin_data[t,]) / fin_data[t,]
    }
    ret_frame <- ret_frame[1:tot,]
  }
  else if (type == 3) {
    for (t in 1:tot) {
      ret_frame[t,] = log(fin_data[t+1,]) - log(fin_data[t,])
    }
    ret_frame <- ret_frame[1:tot,]
  }
  else {
    stop("Invalid type of return")
  }
  if (frame == TRUE) {
    return (ret_frame)
  }
  rLOL <- list()
  for (t in 1:tot) {
    time_i = c();
    for (n in 1:stocks) {
      time_i = c(time_i, ret_frame[t,n])
    }
    rLOL[[t]] <- time_i
  }
  return(rLOL)
}

#---------------------------------------------------------------------------------------------------------------------

#' FFGARCH Model Fitting
#' 
#' @param stocks Matrix or data frame of stock prices (T-times X N-Stocks)
#' @param type Specification of type of stock return (1 for abs, 2 for rel, 3 for log)
#' @param sig_init Initial sigma^2 values for N stocks \code c(sig_11^2, sig_12^2,...,sig_1n^2)
#' @param iter Number of iterations per chain in stan call
#' @param chain Number of chains in stan call
#' @param control_p Optional control par. specs for sampler - list(adapt_delta, max_treedepth)

FFGARCH.fit <- function(stocks, type, sig_init, iter, chain, control_p){
  expect_that(stocks,is_a("data.frame"))
  expect_that(type,is_a("numeric"))
  expect_that(sig_init,is_a("numeric"))
  expect_that(iter,is_a("numeric"))
  expect_that(chain,is_a("numeric"))
  expect_equal(type > 0, TRUE)
  expect_equal(type < 4, TRUE)
  expect_equal(iter > 0, TRUE)
  expect_equal(chain > 0, TRUE)
  formated <- returns_prices(stocks, type, FALSE)
  t <- length(formated)
  n <- length(formated[[1]])
  stan_input <- list (T = t, N = n, ret = formated, sigma_init = sig_init)
  if (!missing(control_p)) {
    FFGARCH_fit <- stan(file = 'FFGARCH.stan', data = stan_input, iter = iter, chains = chain, control = control_p) 
  }
  else {
    FFGARCH_fit <- stan(file = 'FFGARCH.stan', data = stan_input, iter = iter, chains = chain) 
  }
  return(FFGARCH_fit)
}
 
#---------------------------------------------------------------------------------------------------------------------

#' Parameter extraction from FFGARCH object
#' 
#'  @param ffobj FFGARCH-model object

FFGARCH.par <- function(ffobj) {
  expect_that(ffobj, is_a("stanfit"))
  
  Pars <- extract(ffobj)
  data_par <- list(mu = apply(Pars$mu, 2, mean),
                   alpha = apply(Pars$alpha, 2, mean),
                   beta = apply(Pars$beta, 2, mean),
                   gamma = apply(Pars$gamma, 2, mean),
                   W = apply(Pars$W, c(2,3), mean))
  return(data_par)
}

#---------------------------------------------------------------------------------------------------------------------
#' Residuals for all assets at all time intervals
#' 
#' @param ffobj FFGARCH object created by FFGARCH.fit()
#' @param stocks Initial data frame of stock prices
#' @param type Type of financial returns called by FFGARCH.fit() (1 = absolute, 2 = relative, 3 = log)

FFGARCH.resids <- function(ffobj, stocks, type) {
  expect_that(ffobj, is_a("stanfit"))
  expect_that(stocks,is_a("data.frame"))
  expect_that(type,is_a("numeric"))
  expect_equal(type > 0, TRUE)
  expect_equal(type < 4, TRUE)
  
  rets_mod <- apply(extract(ffobj)$ret_out, c(2,3), mean)
  rets_cal <- returns_prices(stocks, type, TRUE)
  return(rets_cal - rets_mod)
}

#---------------------------------------------------------------------------------------------------------------------