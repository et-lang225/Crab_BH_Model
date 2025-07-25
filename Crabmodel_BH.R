C_df <- read.csv('//wlf-statewide.swe.la.gov/FS_WLF/Fisheries/R&A/Age, Growth, & Reproduction/Assessment_Erik/Blue Crab/2022/Model/Crab_Data_BH.csv', header = TRUE)

Bev_Holt <- function(a, b, SSB, recruit_dev = rep(0, length(SSB))){
  SSB_lag <- dplyr::lag(SSB)
  SSB_lag[1] <- SSB[1]
  R <- a*SSB_lag/(b+SSB_lag) * exp(recruit_dev)
  return(R)
}

safe_log <- function(x, min_val = 1e-6) {
  x_clean <- ifelse(is.na(x) | !is.finite(x) | x < min_val, min_val, x)
  return(log(x_clean))
}

S_R_Model <- function(a, b, M, q_r, q_n, catch, r_obs, n_obs, CV_catch, CV_r, CV_n, prop_mat, recruit_dev){
  #Translate CVs into standard errors?
  SE_catch <- safe_log(1+CV_catch^2)
  SE_r <- safe_log(1+CV_r^2)
  SE_n <- safe_log(1+CV_n^2)
  SE_r <- pmax(SE_r, 1e-6)
  SE_n <- pmax(SE_n, 1e-6)
  SE_catch <- pmax(SE_catch, 1e-6)
  # define SSB and R population scale values
  N_obs <- n_obs/(q_n*exp(-0.5*M))
  SSB <- N_obs*prop_mat
  R_obs <- r_obs/(q_r*exp(-0.5*M))
  # Calculate expected values from Beverton Holt
  R_exp <- Bev_Holt(a,b,SSB, recruit_dev)
  R_sq_er <- (log(R_obs) - log(R_exp))^2
  N_exp <- rep(NA, length(N_obs))
  N_exp[1] <- mean(N_obs[1:3], na.rm=TRUE)
  for (i in 2:length(N_obs)){
    if (i==2){
      N_exp[i] = N_exp[i-1]*exp(-M)+R_obs[i-1]*exp(-M)
    }
    else {
      N_exp[i] = N_exp[i-1]*exp(-M)+R_exp[i-1]*exp(-M)
    }
  }
  r_exp <- q_r*R_exp*exp(-0.5*M)
  n_exp <- q_n*N_exp*exp(-0.5*M)
  r_sq_er <- (log(r_obs) - log(r_exp))^2/SE_r^2
  n_sq_er <- (log(n_obs) - log(n_exp))^2/SE_n^2
  # Now catch
  N_obs_lead <- dplyr::lead(N_obs)
  N_obs_lead[length(N_obs_lead)] <- N_obs[length(N_obs)]
  Z = safe_log((R_obs+N_obs)/N_obs_lead)
  F_raw <- Z - M
  F <- ifelse(F_raw < 0, 1e-10, F_raw)
  C_exp <- (F / Z) * N_exp * (1 - exp(-Z))
  C_sq_er <- (log(catch) - log(C_exp))^2/SE_catch^2
  result <- data.frame(r_obs, r_exp, r_sq_er, n_obs, n_exp, n_sq_er, catch, C_exp, C_sq_er, R_sq_er, Z, F, SSB, R_obs, N_obs, R_exp, recruit_dev, N_exp)
  return(result)
}

h <- 0.99
M <- 1

obj_fn <- function(par){
  Rv <- par[1]
  SSBv <- par[2]
  log_q_r <- par[3]
  log_q_n <- par[4]
  recruit_dev <- par[5:length(par)]
  
  result <- S_R_Model(a = (4*h*Rv)/(5*h-1), 
                     b = (SSBv*(1-h))/(5*h-1), 
                     M, 
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
  RDev_penalty <-  sum(recruit_dev^2)
  
  cat("\nn:", n_term, " r:", r_term, " c:", c_term, " R:", R_term, " pen:", RDev_penalty, "\n")
  
  total <- n_term + r_term + c_term + R_term + RDev_penalty
  if (!is.finite(total)) stop("Non-finite obj_fn value detected.")
  return(total)
  }

n_years <- nrow(C_df)
start_dev <- log(C_df$r_obs/mean(C_df$r_obs))
start_par <- c(Rv=600, SSBv=100, log_q_r=-5, log_q_n=-5, R_Devs=start_dev)
lower_bounds <- c(400, 90, -12, -12, rep(-1, n_years))
upper_bounds <- c(800, 200, -2, -2, rep(5, n_years))

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


S_R_Model_fit <- S_R_Model(a_fit,b_fit,M,qr_fit,qn_fit, C_df$Catch_total, C_df$r_obs, C_df$n_obs, C_df$Catch_CV, C_df$CV_r, C_df$CV_n, C_df$adults, recruit_dev_fit)
S_R_Model_fit$SPR <- S_R_Model_fit$SSB/SSBv_fit
SSB_lag <- dplyr::lag(S_R_Model_fit$SSB)

library(ggplot2)
ggplot()+
  ylim(0,max(S_R_Model_fit$R_obs))+
  geom_line(aes(x=S_R_Model_fit$SSB, y=a_fit*SSB_lag/(b_fit+SSB_lag)))+
  geom_point(aes(x=S_R_Model_fit$SSB, y=dplyr::lead(S_R_Model_fit$R_obs)))+
  labs(x='Spawning Stock Biomass', y='Recruitment (1000s of crabs)')+
  theme_bw()

ggplot()+
  ylim(0, max(S_R_Model_fit$SPR))+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$SPR))+
  geom_hline(aes(yintercept=0.3), color='red')+
  labs(x='Year', y='Spawning Potential Ratio')+
  theme_bw()

ggplot()+
  ylim(0, max(C_df$r_obs))+
  geom_point(aes(x=C_df$Year, y=C_df$r_obs))+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$r_exp))+
  theme_bw()

ggplot()+
  ylim(0, max(S_R_Model_fit$n_exp))+
  geom_point(aes(x=C_df$Year, y=C_df$n_obs))+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$n_exp))+
  theme_bw()

ggplot()+
  ylim(0, max(C_df$Catch_total, S_R_Model_fit$C_exp))+
  geom_point(aes(x=C_df$Year, y=C_df$Catch_total), color='blue')+
  geom_line(aes(x=C_df$Year, y=S_R_Model_fit$C_exp))+
  theme_bw()
