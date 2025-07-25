C_df <- read.csv('Crab_Data_BH.csv', header = TRUE)

Bev_Holt <- function(a, b, SSB, recruit_dev = rep(0, length(SSB))){
  SSB_lag <- dplyr::lag(SSB)
  SSB_lag[1] <- SSB[1]
  R <- a*SSB_lag/(b+SSB_lag) * exp(recruit_dev)
  return(R)
}

S_R_Model <- function(h, SSBv, Rv, q_r, q_n, catch, r_obs, n_obs, CV_catch, CV_r, CV_n, prop_mat, recruit_dev){
  #Translate CVs into standard errors?
  SE_catch <- sqrt(log(1+CV_catch^2))
  SE_r <- sqrt(log(1+CV_r^2))
  SE_n <- sqrt(log(1+CV_n^2))
  SE_r <- pmax(SE_r, 1e-6)
  SE_n <- pmax(SE_n, 1e-6)
  SE_catch <- pmax(SE_catch, 1e-6)
  # define SSB and R population scale values
  N_obs <- n_obs/q_n
  SSB <- N_obs*prop_mat
  R_obs <- r_obs/q_r
  N_lead <- dplyr::lead(N_obs)
  N_lead[length(N_lead)] <- N_obs[length(N_obs)]
  Z <- log((R_obs+N_obs)/N_lead)
  u <- catch/(R_obs+N_obs)
  AM <- 1-exp(-Z)
  F <- (u*Z)/AM
  M <- Z-F
  # Calculate expected values from Beverton Holt
  a <- (4*h*Rv)/(5*h-1)
  b <- (SSBv*(1-h))/(5*h-1)
  R_exp <- Bev_Holt(a,b,SSB, recruit_dev)
  R_exp[1] <- mean(R_obs[1:3], na.rm=TRUE)
  R_sq_er <- (log(R_obs) - log(R_exp))^2 / SE_r^2
  N_exp <- rep(NA, length(N_obs))
  N_exp[1] <- mean(N_obs[1:3], na.rm=TRUE)
  for (i in 2:length(N_exp)){
      N_exp[i] = N_exp[i-1]*exp(-Z[i])+R_exp[i-1]*exp(-Z[i])
  }
  r_exp <- q_r*R_exp
  n_exp <- q_n*N_exp
  r_sq_er <- (log(r_obs) - log(r_exp))^2 / SE_r^2
  n_sq_er <- (log(n_obs) - log(n_exp))^2 / SE_n^2
  # Now catch
  C_exp <- (F / Z) * N_exp * (1 - exp(-Z))
  C_sq_er <- (log(catch) - log(C_exp))^2 / SE_catch^2
  result <- data.frame(r_obs, r_exp, r_sq_er, n_obs, n_exp, n_sq_er, catch, C_exp, C_sq_er, R_sq_er, Z, F, M, SSB, R_obs, N_obs, R_exp, N_exp, recruit_dev)
  return(result)
}

h <- 0.99

# Natural Mortality
#Exploitable_Age <- 3 # years
#Hoenig <- exp(1.44-0.982*log(Exploitable_Age))
#Then <- 4.899*Exploitable_Age^-0.916
#Hamel <- 5.4/Exploitable_Age
#M <- 0.5

obj_fn <- function(par){
  Rv <- par[1]
  SSBv <- par[2]
  log_q_r <- par[3]
  log_q_n <- par[4]
  recruit_dev <- par[5:length(par)]
  
  result <- S_R_Model(h=h,
                     Rv=Rv,
                     SSBv=SSBv,  
                     q_r = exp(log_q_r),
                     q_n = exp(log_q_n),
                     catch = C_df$Catch_total,
                     r_obs = C_df$r_obs,
                     n_obs = C_df$n_obs,
                     CV_catch = C_df$Catch_CV,
                     CV_r = C_df$CV_r,
                     CV_n = C_df$CV_n,
                     prop_mat = C_df$adults,
                     recruit_dev = recruit_dev)
  
  #Total objective to minimize sum of all error terms
  n_term <- sum(result$n_sq_er)
  r_term <- sum(result$r_sq_er)
  c_term <- sum(result$C_sq_er)
  R_term <- sum(result$R_sq_er)
  delta_dev <- diff(recruit_dev)
  dev_penalty <- sum(delta_dev^2, na.rm=TRUE)
  w_dev <- 1  # try different weights
  
  cat("\nn:", n_term, " r:", r_term, " c:", c_term, " R:", R_term, "pen:", dev_penalty)
  
  w_n <- 1
  w_r <- 1
  w_c <- 1
  w_BH <- 1

  total <- w_n * n_term + w_r * r_term + w_c * c_term + w_BH * R_term + w_dev*dev_penalty
  if (!is.finite(total)) stop("Non-finite obj_fn value detected.")
  return(total)
  }

n_years <- nrow(C_df)
start_dev <- log(C_df$r_obs/mean(C_df$r_obs))
start_par <- c(Rv=600, SSBv=100, log_q_r=-5, log_q_n=-5, R_Devs=start_dev)
lower_bounds <- c(400, 90, -10, -10, rep(-10, n_years))
upper_bounds <- c(700, 200, -2, -2, rep(10, n_years))

fit <- optim(start_par, obj_fn, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, control=list(maxit=1000))
fit$convergence
fit$message
fit$value
fit$counts
fit$par

Rv_fit <- fit$par['Rv']
SSBv_fit <- fit$par['SSBv']
qr_fit <- exp(fit$par['log_q_r'])
qn_fit <- exp(fit$par['log_q_n'])
recruit_dev_fit <- fit$par[5:length(fit$par)]

a_fit <- (4*h*Rv_fit)/(5*h-1)
b_fit <- (SSBv_fit*(1-h))/(5*h-1)


S_R_Model_fit <- S_R_Model(h, SSBv_fit, Rv_fit, qr_fit, qn_fit, C_df$Catch_total, C_df$r_obs, C_df$n_obs, C_df$Catch_CV, C_df$CV_r, C_df$CV_n, C_df$adults, recruit_dev_fit)
S_R_Model_fit$SPR <- S_R_Model_fit$SSB/SSBv_fit
S_R_Model_fit$Mort <- exp(-log((S_R_Model_fit$R_exp+S_R_Model_fit$N_exp)/dplyr::lead(S_R_Model_fit$N_exp)))
SSB_lag <- dplyr::lag(S_R_Model_fit$SSB)


library(ggplot2)
ggplot()+
  ylim(0, max(C_df$r_obs, S_R_Model_fit$r_exp))+
  geom_point(aes(x=C_df$Year, y=C_df$r_obs))+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$r_exp))+
  theme_bw()

ggplot()+
  ylim(0, max(C_df$n_obs,S_R_Model_fit$n_exp))+
  geom_point(aes(x=C_df$Year, y=C_df$n_obs))+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$n_exp))+
  theme_bw()

ggplot()+
  ylim(0, max(C_df$Catch_total, S_R_Model_fit$C_exp))+
  geom_point(aes(x=C_df$Year, y=C_df$Catch_total), color='blue')+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$C_exp))+
  theme_bw()

ggplot()+
geom_point(aes(x=C_df$Year, y=S_R_Model_fit$recruit_dev))+
geom_hline(aes(yintercept=0), color='red')+
theme_bw()

SSB_series <- seq(0,max(S_R_Model_fit$SSB), by=10)
BH_mean <- a_fit*SSB_series/(b_fit+SSB_series)
ggplot()+
  ylim(0,max(S_R_Model_fit$R_obs, BH_mean))+
  geom_line(aes(x=SSB_series, y=BH_mean))+
  geom_point(aes(x=S_R_Model_fit$SSB, y=dplyr::lead(S_R_Model_fit$R_obs)))+
  labs(x='Spawning Stock Biomass', y='Recruitment (1000s of crabs)')+
  theme_bw()

ggplot()+
  ylim(0, max(S_R_Model_fit$SPR))+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$SPR))+
  geom_hline(aes(yintercept=0.3), color='red')+
  labs(x='Year', y='Spawning Potential Ratio')+
  theme_bw()

