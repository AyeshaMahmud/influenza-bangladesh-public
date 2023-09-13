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
library(HDInterval)
require(lhs)


#### set parameters for tuning ####
rm(list = ls())

#### round 1 of tuning ####

parms <- readRDS(here("data", "cleaned", "LHS_5000_parameter_matrix_NULL2.RDS"))

districts_list_limited <- readRDS(here("data", "cleaned", "district_list_limited.RDS"))

# select top 1000 parameter combinations based on RMSE
tuning_round = 1
perc_change_bound = 10
num_ens_keep = 500
num_sim = 1
parms_df <- parms[1, ] %>% mutate(District = NA)



for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_NULL2.RDS"))) #original simulation
  mean_rmse_df <- tmp[[3]] 
  top_ens <- mean_rmse_df$num_ens[1:num_ens_keep]
  

  parms_tmp <- parms[top_ens, ] %>% mutate(District = district_name)
  parms_df <- rbind(parms_df, parms_tmp)

}

parms_df <- parms_df %>% filter(!is.na(District))      
parms_df_long <- parms_df %>% gather(parameter, value, -c(District))

# for each parameter find the 95% HDI and create new list of parameter ranges

hdi_range <- as.data.frame(cbind(lower = 1, upper = 1, parameter = NA))

varlist =  c("R0max","lambda","D", "L", "S0perc", "I0")
for (v in 1:length(varlist)){
 lower <- hdi(parms_df[,eval(varlist[v])])[1] %>% unlist() %>% as.numeric() %>% round(.,2)
 upper <- hdi(parms_df[,eval(varlist[v])])[2] %>% unlist() %>% as.numeric() %>% round(.,2)
 hdi_range <- rbind.data.frame(hdi_range, cbind.data.frame(lower = as.numeric(lower), upper = as.numeric(upper), parameter = varlist[v]))
}

hdi_range <- hdi_range %>% filter(!is.na(parameter))
lowerbound <- c(1.5, 0.1, 2, 365, 0.4, 500) #original bounds
upperbound <- c(3, 0.8, 5, 356*4, 0.8, 1000)
names(lowerbound) <- names(upperbound) <- c("R0max", "lambda","D", "L", "S0perc", "I0")

hdi_range$lower_perc_change <- round((hdi_range$lower - lowerbound)*100/lowerbound) 
hdi_range$upper_perc_change <- round((upperbound - hdi_range$upper)*100/upperbound)
hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound][!is.na(hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound])] #D
hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound] #None

saveRDS(hdi_range, file = here("data", "cleaned", paste0("hdi_range_tuning_round",tuning_round,"_NULL2.RDS")))


# update lower and upper bound
lowerbound[hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound]] <- hdi_range$lower[hdi_range$lower_perc_change >= perc_change_bound]
upperbound[hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound]]<- hdi_range$upper[hdi_range$upper_perc_change >= perc_change_bound]
lowerbound <- lowerbound[!is.na(lowerbound)]
upperbound <- upperbound[!is.na(upperbound)]

saveRDS(list(lowerbound, upperbound), file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round,"_parameter_matrix_updatedbounds_NULL2.RDS")))


#### generate paramter combinations using new ranges
num_ens <- 5000
set.seed(123456)
x <- randomLHS(num_ens, 6) 
y <- x 

y[,1] <- qunif(x[,1], 1.5, 3) #R0max 
y[,2] <- qunif(x[,2], 0.1, 0.8) #lambda #amplitude
y[,3] <- qunif(x[,3], hdi_range$lower[hdi_range$parameter == "D"], 5) #D
y[,4] <- qunif(x[,4], 365, 365*4) #L #wider interval used here
y[,5] <- qunif(x[,5], 0.4, 0.8) #S0perc
y[,6] <- qunif(x[,6], 500, 1000) #I0

parms <- y %>% as.data.frame()
names(parms) <- c("R0max", "lambda", "D", "L", "S0perc", "I0")
saveRDS(parms, file = here("data", "cleaned", paste0("LHS_",num_ens,"_parameter_matrix_tuning_round",tuning_round,"_NULL2.RDS")))


#### round 2 of tuning ####
tuning_round = 2
perc_change_bound = 10
num_ens_keep = 500
parms <- readRDS(here("data", "cleaned", paste0("LHS_5000_parameter_matrix_tuning_round",tuning_round - 1,"_NULL2.RDS")))
districts_list_limited <- readRDS(here("data", "cleaned", "district_list_limited.RDS"))

# select top 1000 parameter combinations based on RMSE
num_sim = 1
parms_df <- parms[1, ] %>% mutate(District = NA)


for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_tuning_round",tuning_round,"_NULL2.RDS")))
  mean_rmse_df <- tmp[[3]] 
  top_ens <- mean_rmse_df$num_ens[1:num_ens_keep]
  
  
  parms_tmp <- parms[top_ens, ] %>% mutate(District = district_name)
  parms_df <- rbind(parms_df, parms_tmp)
  
}

parms_df <- parms_df %>% filter(!is.na(District))      
parms_df_long <- parms_df %>% gather(parameter, value, -c(District))

# for each parameter find the 95% HDI and create new list of parameter ranges
hdi_range <- as.data.frame(cbind(lower = 1, upper = 1, parameter = NA))

varlist =  c("R0max","lambda","D", "L", "S0perc", "I0")
for (v in 1:length(varlist)){
  lower <- hdi(parms_df[,eval(varlist[v])])[1] %>% unlist() %>% as.numeric() %>% round(.,2)
  upper <- hdi(parms_df[,eval(varlist[v])])[2] %>% unlist() %>% as.numeric() %>% round(.,2)
  hdi_range <- rbind.data.frame(hdi_range, cbind.data.frame(lower = as.numeric(lower), upper = as.numeric(upper), parameter = varlist[v]))
}

hdi_range <- hdi_range %>% filter(!is.na(parameter))
lowerbound = readRDS(file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round-1,"_parameter_matrix_updatedbounds_NULL2.RDS")))[[1]]
upperbound = readRDS(file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round-1,"_parameter_matrix_updatedbounds_NULL2.RDS")))[[2]]

hdi_range$lower_perc_change <- round((hdi_range$lower - lowerbound)*100/lowerbound) 
hdi_range$upper_perc_change <- round((upperbound - hdi_range$upper)*100/upperbound)
hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound][!is.na(hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound])] #R0max, L
hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound] #None

saveRDS(hdi_range, file = here("data", "cleaned", paste0("hdi_range_tuning_round",tuning_round,"_NULL2.RDS")))

# update lower and upper bound
lowerbound[hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound]] <- hdi_range$lower[hdi_range$lower_perc_change >= perc_change_bound]
upperbound[hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound]]<- hdi_range$upper[hdi_range$upper_perc_change >= perc_change_bound]
lowerbound <- lowerbound[!is.na(lowerbound)]
upperbound <- upperbound[!is.na(upperbound)]

saveRDS(list(lowerbound, upperbound), file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round,"_parameter_matrix_updatedbounds_NULL2.RDS")))


#### generate paramter combinations using new ranges
num_ens <- 5000 
x <- randomLHS(num_ens, 6) 
y <- x 

y[,1] <- qunif(x[,1], hdi_range$upper[hdi_range$parameter == "R0max"], 3) #R0max 
y[,2] <- qunif(x[,2], 0.1, 0.8) #lambda #amplitude
y[,3] <- qunif(x[,3], hdi_range$lower[hdi_range$parameter == "D"], 5) #D
y[,4] <- qunif(x[,4], hdi_range$upper[hdi_range$parameter == "L"], 365*4) #L 
y[,5] <- qunif(x[,5], 0.4, 0.8) #S0perc
y[,6] <- qunif(x[,6], 500, 1000) #I0

parms <- y %>% as.data.frame()
names(parms) <- c("R0max", "lambda", "D", "L", "S0perc", "I0")
saveRDS(parms, file = here("data", "cleaned", paste0("LHS_",num_ens,"_parameter_matrix_tuning_round",tuning_round,"_NULL2.RDS")))


#### round 3 of tuning ####
tuning_round = 3
perc_change_bound = 10
num_ens_keep = 500
parms <- readRDS(here("data", "cleaned", paste0("LHS_5000_parameter_matrix_tuning_round",tuning_round - 1,"_NULL2.RDS")))
districts_list_limited <- readRDS(here("data", "cleaned", "district_list_limited.RDS"))

# select top 1000 parameter combinations based on RMSE
num_sim = 1
parms_df <- parms[1, ] %>% mutate(District = NA)


for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_tuning_round",tuning_round,"_NULL2.RDS")))
  mean_rmse_df <- tmp[[3]] 
  top_ens <- mean_rmse_df$num_ens[1:num_ens_keep]
  
  
  parms_tmp <- parms[top_ens, ] %>% mutate(District = district_name)
  parms_df <- rbind(parms_df, parms_tmp)
  
}

parms_df <- parms_df %>% filter(!is.na(District))      
parms_df_long <- parms_df %>% gather(parameter, value, -c(District))

# for each parameter find the 95% HDI and create new list of parameter ranges
hdi_range <- as.data.frame(cbind(lower = 1, upper = 1, parameter = NA))

varlist =  c("R0max","lambda","D", "L", "S0perc", "I0")
for (v in 1:length(varlist)){
  lower <- hdi(parms_df[,eval(varlist[v])])[1] %>% unlist() %>% as.numeric() %>% round(.,2)
  upper <- hdi(parms_df[,eval(varlist[v])])[2] %>% unlist() %>% as.numeric() %>% round(.,2)
  hdi_range <- rbind.data.frame(hdi_range, cbind.data.frame(lower = as.numeric(lower), upper = as.numeric(upper), parameter = varlist[v]))
}

hdi_range <- hdi_range %>% filter(!is.na(parameter))
lowerbound = readRDS(file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round-1,"_parameter_matrix_updatedbounds_NULL2.RDS")))[[1]]
upperbound = readRDS(file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round-1,"_parameter_matrix_updatedbounds_NULL2.RDS")))[[2]]

hdi_range$lower_perc_change <- round((hdi_range$lower - lowerbound)*100/lowerbound) 
hdi_range$upper_perc_change <- round((upperbound - hdi_range$upper)*100/upperbound)
hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound][!is.na(hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound])] #R0max, L
hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound] #None

saveRDS(hdi_range, file = here("data", "cleaned", paste0("hdi_range_tuning_round",tuning_round,"_NULL2.RDS")))

# update lower and upper bound
lowerbound[hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound]] <- hdi_range$lower[hdi_range$lower_perc_change >= perc_change_bound]
upperbound[hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound]]<- hdi_range$upper[hdi_range$upper_perc_change >= perc_change_bound]
lowerbound <- lowerbound[!is.na(lowerbound)]
upperbound <- upperbound[!is.na(upperbound)]

saveRDS(list(lowerbound, upperbound), file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round,"_parameter_matrix_updatedbounds_NULL2.RDS")))


#### generate paramter combinations using new ranges
num_ens <- 5000 
x <- randomLHS(num_ens, 6) 
y <- x 

y[,1] <- qunif(x[,1], hdi_range$upper[hdi_range$parameter == "R0max"], 3) #R0max 
y[,2] <- qunif(x[,2], 0.1, 0.8) #lambda #amplitude
y[,3] <- qunif(x[,3], hdi_range$lower[hdi_range$parameter == "D"], 5) #D
y[,4] <- qunif(x[,4], hdi_range$upper[hdi_range$parameter == "L"], 365*4) #L 
y[,5] <- qunif(x[,5], 0.4, 0.8) #S0perc
y[,6] <- qunif(x[,6], 500, 1000) #I0

parms <- y %>% as.data.frame()
names(parms) <- c("R0max", "lambda", "D", "L", "S0perc", "I0")
saveRDS(parms, file = here("data", "cleaned", paste0("LHS_",num_ens,"_parameter_matrix_tuning_round",tuning_round,"_NULL2.RDS")))

#### round 4 of tuning ####
tuning_round = 4
perc_change_bound = 10
num_ens_keep = 500
parms <- readRDS(here("data", "cleaned", paste0("LHS_5000_parameter_matrix_tuning_round",tuning_round - 1,"_NULL2.RDS")))
districts_list_limited <- readRDS(here("data", "cleaned", "district_list_limited.RDS"))

# select top 1000 parameter combinations based on RMSE
num_sim = 1
parms_df <- parms[1, ] %>% mutate(District = NA)


for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_tuning_round",tuning_round,"_NULL2.RDS")))
  mean_rmse_df <- tmp[[3]] 
  top_ens <- mean_rmse_df$num_ens[1:num_ens_keep]
  
  
  parms_tmp <- parms[top_ens, ] %>% mutate(District = district_name)
  parms_df <- rbind(parms_df, parms_tmp)
  
}

parms_df <- parms_df %>% filter(!is.na(District))      
parms_df_long <- parms_df %>% gather(parameter, value, -c(District))

# for each parameter find the 95% HDI and create new list of parameter ranges
hdi_range <- as.data.frame(cbind(lower = 1, upper = 1, parameter = NA))

varlist =  c("R0max","lambda","D", "L", "S0perc", "I0")
for (v in 1:length(varlist)){
  lower <- hdi(parms_df[,eval(varlist[v])])[1] %>% unlist() %>% as.numeric() %>% round(.,2)
  upper <- hdi(parms_df[,eval(varlist[v])])[2] %>% unlist() %>% as.numeric() %>% round(.,2)
  hdi_range <- rbind.data.frame(hdi_range, cbind.data.frame(lower = as.numeric(lower), upper = as.numeric(upper), parameter = varlist[v]))
}

hdi_range <- hdi_range %>% filter(!is.na(parameter))
lowerbound = readRDS(file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round-1,"_parameter_matrix_updatedbounds_NULL2.RDS")))[[1]]
upperbound = readRDS(file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round-1,"_parameter_matrix_updatedbounds_NULL2.RDS")))[[2]]

hdi_range$lower_perc_change <- round((hdi_range$lower - lowerbound)*100/lowerbound) 
hdi_range$upper_perc_change <- round((upperbound - hdi_range$upper)*100/upperbound)
hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound][!is.na(hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound])] #R0max, L
hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound] #None

saveRDS(hdi_range, file = here("data", "cleaned", paste0("hdi_range_tuning_round",tuning_round,"_NULL2.RDS")))

# update lower and upper bound
lowerbound[hdi_range$parameter[hdi_range$lower_perc_change >= perc_change_bound]] <- hdi_range$lower[hdi_range$lower_perc_change >= perc_change_bound]
upperbound[hdi_range$parameter[hdi_range$upper_perc_change >= perc_change_bound]]<- hdi_range$upper[hdi_range$upper_perc_change >= perc_change_bound]
lowerbound <- lowerbound[!is.na(lowerbound)]
upperbound <- upperbound[!is.na(upperbound)]

saveRDS(list(lowerbound, upperbound), file = here("data", "cleaned", paste0("LHS_tuning_round",tuning_round,"_parameter_matrix_updatedbounds_NULL2.RDS")))


#### generate paramter combinations using new ranges
num_ens <- 5000*2 
x <- randomLHS(num_ens, 6) 
y <- x 

y[,1] <- qunif(x[,1], hdi_range$upper[hdi_range$parameter == "R0max"], 3) #R0max 
y[,2] <- qunif(x[,2], 0.1, 0.8) #lambda #amplitude
y[,3] <- qunif(x[,3], 2, 5) #D
y[,4] <- qunif(x[,4], hdi_range$upper[hdi_range$parameter == "L"], 365*4) #L 
y[,5] <- qunif(x[,5], 0.4, 0.8) #S0perc
y[,6] <- qunif(x[,6], 500, 1000) #I0

parms <- y %>% as.data.frame()
names(parms) <- c("R0max", "lambda", "D", "L", "S0perc", "I0")
saveRDS(parms, file = here("data", "cleaned", paste0("LHS_",num_ens,"_parameter_matrix_tuning_round",tuning_round,"_NULL2.RDS")))

