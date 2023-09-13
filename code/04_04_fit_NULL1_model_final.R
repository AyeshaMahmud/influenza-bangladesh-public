library(lubridate)
library(stringr)
library(tidyverse)
library(dplyr)
library(here)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(here)
library(sjPlot)
library(viridis)
library(cowplot)
library(sandwich)
library(stringr)
library(gridExtra)
library(grid)
library(lattice)
library(ggpubr)
library(wesanderson)
library(scales)


rm(list = ls())

# this is the final run of the NULL1 model using 10,000 paramter combinations using bounds after parameter tuning

#### set up ####
source(here('code', 'AH_T_model.R'))

districts_list <- readRDS(here("data", "cleaned", "district_list_limited.RDS")) #districts must have at least five years of data prior to testing period (2018 and 2019)
num_sim = 1

pop = read.csv(file = here("data", "district_pop.csv"), header = TRUE)

# fixed paramters
mu = (18/1000)/365
expoI = 0.97
#N = 100000

#### fit model for each district ####
for (d in 9:length(districts_list)) {
#for (d in c(1,4,5,12)) {
    
  # set name of district
  district_name = districts_list[d]
  
  # set population
  N = pop$N[pop$District == district_name]
  #### get climate data ####
  
  daily_clim <- readRDS(here("data", "cleaned", "era5_clim_data_for_model.RDS")) %>%
    filter(District == district_name)
  
  # create a historical timeseries to deal with transients
  burnin <- 20
  q_ts_historic <- daily_clim  %>% arrange(doy) %>%
    group_by(doy)  %>% dplyr::summarize(SH = mean(as.numeric(SH),na.rm=TRUE)) %>%
    select(SH) %>% unlist() %>% as.numeric() %>% 
    rep(., times = burnin) 
  q_ts <- c(q_ts_historic, daily_clim$SH)
  
  T_ts_historic <- daily_clim %>% 
    mutate(doy = strftime(date, format = "%j")) %>%
    group_by(doy) %>% arrange(doy) %>% dplyr::summarize(tmp = mean(as.numeric(tmp),na.rm=T)) %>%
    select(tmp) %>% unlist() %>% as.numeric() %>% 
    rep(., times = burnin) 
  
  T_ts <- c(T_ts_historic, daily_clim$tmp)
  
  #### observed data ####
   
  obs <- readRDS(here("data", "cleaned", "observed_monthly_dat_for_model.RDS")) %>% filter(District == district_name)  %>% arrange(date)
  
  #### parameter matrix ####
  parms <- readRDS(here("data", "cleaned", paste0("LHS_10000_parameter_matrix_tuning_round",4,"_NULL1.RDS")))
  parms <- parms %>% select('R0max','D','L','S0perc','I0')
  num_ens=dim(parms)[1]
  R0 = matrix(parms[,1], length(q_ts), num_ens, byrow = TRUE) # set at fixed value
  
  
  # Get other params:
  D <- parms[, 'D']; S0perc <- parms[, 'S0perc']; I0 <- parms[, 'I0']; L <- parms[,"L"]
  
  
  S0 = S0perc*N
  tm_strt <- 1; tm_end <- length(q_ts) 
  tm_step <- 1; tm.range <- tm_strt:tm_end
  
  # Calculate betas:
  beta <- sapply(1:num_ens, function(ix) {
    R0[, ix] / D[ix]
  })
  rm(R0)
  
  
  
  
  date_vec <- daily_clim$date
  
  if(num_sim > 1) {discrete = TRUE} else {discrete = FALSE}
  
  rmse <- matrix(NA, num_ens, num_sim)
  
  #### simulate dynamics for num_sim number of rounds #### 
  #if stochastic, num_sim should be greater than 1
  for (s in 1:num_sim) {
  
    if(discrete == TRUE) {set.seed(s)} #set seed so we can replicate
    sim_AH_T <- SIRS_AH_T(tm_strt, tm_end, tm_step, S0, I0, N, D, L, mu, beta, expoI, realdata=TRUE)
    inc <- sapply(1:num_ens, function(ix) {
      diff(c(0,sim_AH_T$newI[,ix]))[1:length(sim_AH_T$newI[,ix])]
    })
    
    
   out <- inc[(length(q_ts_historic)+2):length(sim_AH_T$newI[,1]), ] %>% as.data.frame() %>% gather(num_ens, inc) %>% mutate(num_ens = parse_number(num_ens))
   out$date <- rep(date_vec, length.out = nrow(out))
   out$month <- month(out$date)
   out$year <- year(out$date)
   
   out_monthly_all <- out  %>%
     group_by(num_ens, year, month) %>% summarize(total_inc = sum(inc, na.rm = T), .groups = 'drop') %>%
     left_join(., obs, by = c("year", "month")) %>%
     ungroup() %>% 
     mutate(total_inc_missing = case_when(year <= 2017 & !is.na(mean) ~ as.numeric(total_inc),
                                          year <= 2017 & is.na(mean) ~ mean,
                                          year > 2017 ~ as.numeric(NA))) %>%
     group_by(num_ens) %>%
     mutate(sim_scaled = rescale(total_inc_missing)) 
   
   # fit data from 2012 to 2017
   rmse_df <- out_monthly_all %>% filter(year <= 2017) %>% mutate(error = (obs_monthly_scaled - sim_scaled)^2) %>%
     group_by(num_ens) %>% 
     summarize(sse = sum(error, na.rm = TRUE), n = length(!is.na(obs_monthly_scaled)),  .groups = 'drop') %>%
     mutate(rmse = sqrt(sse/n))
   
   rmse[,s] <- rmse_df$rmse

      set.seed(NULL)
  }
  
  
  mean_rmse <- rowMeans(rmse, na.rm = TRUE)
  mean_rmse_df <- as.data.frame(mean_rmse)
  mean_rmse_df$num_ens <- seq(1:nrow(mean_rmse_df))
  mean_rmse_df <- mean_rmse_df %>% arrange(mean_rmse)
  top_ten_ens <- mean_rmse_df$num_ens[1:10]
  
  
  
  out_df_top_ten <- out_monthly_all %>% filter(num_ens %in% top_ten_ens)
  
  
  saveRDS(list(out_monthly_all, out_df_top_ten, mean_rmse_df, top_ten_ens, out), file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_final_NULL1.RDS")))
  

  
}



