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
library(sf)
library(nngeo)
library(scales)
require(lhs)

rm(list = ls())

## Notes:
# Assuming NISB and HBIS data are extracted and cleaned and saved in data > cleaned folder
# Assuming ERA5 climate data have been downloaded and saved in data > cleaned folder
# Assuming CMIP6 projection data have been downloaded and saved in data > cleaned folder


#### process observed data ####
nisb1 <- readRDS(here("data", "cleaned", list.files(here("data", "cleaned"), pattern="cleaned_NISB_latest_"))[1])
nisb2 <- readRDS(here("data", "cleaned", list.files(here("data", "cleaned"), pattern="cleaned_NISB_latest_"))[2])
hbis1 <- readRDS(here("data","cleaned", list.files(here("data", "cleaned") , pattern="cleaned_HBIS_latest_"))[1])
hbis2 <- readRDS(here("data","cleaned", list.files(here("data", "cleaned") , pattern="cleaned_HBIS_latest_"))[2])

nisb1 <- nisb1 %>% rename(Hospital = Site)
nisb2 <- nisb2 %>% rename(Hospital = Site)

# create full dataset by combining the two sources
df_flu_all <- rbind(nisb1, nisb2, hbis1, hbis2) %>% mutate(month = month(date), year = year(date)) %>%
  arrange(Hospital, year, month) %>%
  mutate(lag_Prop = lag(Prop)) %>%
  filter(District != "Total")

df_flu_all$District[grepl("Barishal",df_flu_all$District, fixed = TRUE)]<- "Barisal"
df_flu_all$District[grepl("Potuakhali",df_flu_all$District, fixed = TRUE)]<- "Patuakhali"
df_flu_all$District[grepl("Kishoregonj",df_flu_all$District, fixed = TRUE)]<- "Kishoreganj"

obs <- df_flu_all %>% group_by(District, date, month, year) %>% filter(year < 2020) %>%
  dplyr::summarize(mean = mean(Prop, na.rm = TRUE)) %>% group_by(District) %>% 
  dplyr::mutate(mean_test_data = case_when(year > 2017 ~ as.numeric(NA),
                                           year <= 2017 ~ as.numeric(mean))) %>% #for model fitting
  dplyr::mutate(obs_monthly_scaled = rescale(mean_test_data, na.rm = TRUE)) %>%
  dplyr::mutate(obs_monthly_scaled_full = rescale(mean, na.rm = TRUE)) %>%
  arrange(District, date, month, year)

saveRDS(obs, here("data", "cleaned", "observed_monthly_dat_for_model.RDS"))

#### get names of all district ####
districts_list <- unique(obs$District)
saveRDS(districts_list, file = here("data", "cleaned", "district_list.RDS"))

#### get names of districts that have at least five years of data for fitting ####
districts_list_limited <- obs %>% group_by(District, year) %>% summarize(count = n()) %>%
  mutate(full_year = case_when(count == 12 ~ 1, TRUE~0)) %>% filter(year < 2018) %>%
  group_by(District) %>% summarize(num_full_year = sum(full_year))  %>% filter(num_full_year >= 5) %>% select(District) %>% unlist()

saveRDS(districts_list_limited, file = here("data", "cleaned", "district_list_limited.RDS"))

#### process ERA5 data on temperature and humidity ####
tmp_dat <- readRDS(here('data', 'cleaned', 'all_tmp_by_site_era5.rds')) %>%
  group_by(District, date, doy, year, month) %>%
  summarize(tmp = mean(tmp, na.rm=TRUE), tmp_K = mean(tmp_K, na.rm = TRUE))


clim_dat <- readRDS(here('data', 'cleaned', 'all_SH_by_site_era5.rds')) %>%
  group_by(District, date, doy, year, month) %>%
  summarize(SH = mean(SH,na.rm = TRUE)) %>%
  left_join(.,tmp_dat, by = c("District", "month", "year", "date", "doy")) 

saveRDS(clim_dat, file = here("data", "cleaned", "era5_clim_data_for_model.RDS"))


#### process CMIP6 projection data on temperature and humidity ####
tmp_dat_proj <- readRDS(here('data', 'cleaned', 'all_tmp_by_site_CMIP6_GFDL.rds')) %>%
  group_by(District, date, doy, year, month) %>%
  summarize(tmp = mean(tmp, na.rm=TRUE), tmp_K = mean(tmp_K, na.rm = TRUE))


clim_dat_proj <- readRDS(here('data', 'cleaned', 'all_SH_by_site_CMIP6_GFDL.rds')) %>%
  group_by(District, date, doy, year, month) %>%
  summarize(SH = mean(SH,na.rm = TRUE)) %>%
  left_join(.,tmp_dat_proj, by = c("District", "month", "year", "date", "doy")) 

saveRDS(clim_dat_proj, file = here("data", "cleaned", "CMIP6_GFDL_clim_data_for_model.RDS"))


# historical CMIP6 data
tmp_dat_proj <- readRDS(here('data', 'cleaned', 'all_tmp_by_site_CMIP6_GFDL_historical.rds')) %>%
  group_by(District, date, doy, year, month) %>%
  summarize(tmp = mean(tmp, na.rm=TRUE), tmp_K = mean(tmp_K, na.rm = TRUE))


clim_dat_proj <- readRDS(here('data', 'cleaned', 'all_SH_by_site_CMIP6_GFDL_historical.rds')) %>%
  group_by(District, date, doy, year, month) %>%
  summarize(SH = mean(SH,na.rm = TRUE)) %>%
  left_join(.,tmp_dat_proj, by = c("District", "month", "year", "date", "doy")) 

saveRDS(clim_dat_proj, file = here("data", "cleaned", "CMIP6_GFDL_clim_data_for_model_historical.RDS"))


# Note: to look at climate change rather than variability, we calculate average by day of year for both baseline and future
# timeseries; take the difference; then add difference to baseline era5 data

# difference between projection and climate by day of year 
baseline_clim <- readRDS(file = here("data", "cleaned", "CMIP6_GFDL_clim_data_for_model_historical.RDS")) %>% 
  filter(year >= 2005 & year <= 2020) %>%
  group_by(doy, District) %>%
  summarize(SH = mean(SH, na.rm = TRUE), tmp = mean(tmp, na.rm = TRUE))

proj_clim <- readRDS(file = here("data", "cleaned", "CMIP6_GFDL_clim_data_for_model.RDS")) %>% 
  filter(year >= 2085 & year <= 2100) %>%
  group_by(doy, District) %>%
  summarize(SH = mean(SH, na.rm = TRUE), tmp = mean(tmp, na.rm = TRUE))


doy_diff <- baseline_clim %>% left_join(.,proj_clim, by = c("doy", "District")) %>%
  mutate(SH_diff = SH.y - SH.x, tmp_diff = tmp.y - tmp.x) %>%
  select(doy, District, SH_diff, tmp_diff)

# create new baseline
new_baseline_clim <- readRDS(file = here("data", "cleaned", "era5_clim_data_for_model.RDS")) %>%
  left_join(., doy_diff, by = c("doy", "District")) %>%
  mutate(SH_cc = SH + SH_diff, tmp_cc = tmp + tmp_diff)

saveRDS(new_baseline_clim, file = here("data", "cleaned", "era5_CMIP6_GFDL_modified_clim_data_for_projection.RDS"))



#### create parameter list using ranges and LHS ####
num_ens <- 5000
set.seed(12345)
x <- randomLHS(num_ens, 11) 
y <- x 

y[,1] <- qunif(x[,1], 1.5, 3) #R0max 
y[,2] <- qunif(x[,2], 0.86, 1.18) #R0diff
y[,3] <- qunif(x[,3], 2.2, 4) #qmin
y[,4] <- qunif(x[,4], 17, 20) #qmax
y[,5] <- qunif(x[,5], 4, 12) #qmid 
y[,6] <- qunif(x[,6], 20.4, 30) #Tc 
#y[,7] <- qunif(x[,7], 0, 15) #Tdiff
y[,7] <- qunif(x[,7], 0.95, 1.54) #Texp
y[,8] <- qunif(x[,8], 2, 5) #D
y[,9] <- qunif(x[,9], 365, 365*4) #L 
y[,10] <- qunif(x[,10], 0.4, 0.8) #S0perc
y[,11] <- qunif(x[,11], 500, 1000) #I0

parms <- y %>% as.data.frame()
names(parms) <- c("R0max", "R0diff", "qmin", "qmax", "qmid", "Tc", "Texp", "D", "L", "S0perc", "I0")
parms$Tdiff <- 0 # not using this parameter
saveRDS(parms, file = here("data", "cleaned", paste0("LHS_",num_ens,"_parameter_matrix.RDS")))


#### create parameter list for NULL 1 model (constant R0) ####
num_ens <- 5000
set.seed(12345)
x <- randomLHS(num_ens, 5) 
y <- x 

y[,1] <- qunif(x[,1], 1.5, 3) #R0max 
y[,2] <- qunif(x[,2], 2, 5) #D
y[,3] <- qunif(x[,3], 365, 365*4) #L 
y[,4] <- qunif(x[,4], 0.4, 0.8) #S0perc
y[,5] <- qunif(x[,5], 500, 1000) #I0

parms <- y %>% as.data.frame()
names(parms) <- c("R0max", "D", "L", "S0perc", "I0")
saveRDS(parms, file = here("data", "cleaned", paste0("LHS_",num_ens,"_parameter_matrix_NULL1.RDS")))



#### create parameter list for NULL 2 model (sinusoidal R0) ####
num_ens <- 5000
set.seed(12345)
x <- randomLHS(num_ens, 6) 
y <- x 

y[,1] <- qunif(x[,1], 1.5, 3) #R0max 
y[,2] <- qunif(x[,2], 0.1, 0.8) #lambda (amplitude)
y[,3] <- qunif(x[,3], 2, 5) #D
y[,4] <- qunif(x[,4], 365, 365*4) #L 
y[,5] <- qunif(x[,5], 0.4, 0.8) #S0perc
y[,6] <- qunif(x[,6], 500, 1000) #I0

parms <- y %>% as.data.frame()
names(parms) <- c("R0max", "lambda","D", "L", "S0perc", "I0")
saveRDS(parms, file = here("data", "cleaned", paste0("LHS_",num_ens,"_parameter_matrix_NULL2.RDS")))



