#### setup ####

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
library(xtable)

## Note: this file produces the final results and figures 

rm(list = ls())

parms <- readRDS(here("data", "cleaned", "LHS_10000_parameter_matrix_tuning_round4.RDS"))
districts_list_limited <- readRDS(here("data", "cleaned", "district_list_limited.RDS"))
source(here('code', 'AH_T_model.R'))
source(here("code", "cg_circular.R"))


#### read in simulation results and store best-fit parameter estimates and simulations by district and get the overall best-fit parameters ####
num_sim = 1
top_ten_df <- as.data.frame(cbind(num_ens = NA, rank = NA, District = NA))
rmse_df <- as.data.frame(cbind(mean_rmse= NA, num_ens = NA, District = NA))


parms_df <- parms[1, ] %>% mutate(District = NA)
sim_df_top_ten<- as.data.frame(cbind(sim_scaled = NA, num_ens = NA, year = NA, month = NA, District = NA, total_inc = NA, mean = NA, mean_test_data = NA, 
                                     obs_monthly_scaled = NA, obs_monthly_scaled_full = NA, total_inc_missing = NA,
                                     total_inc_missing_full = NA, sim_scaled_full = NA)) %>% mutate(date = as.Date("2012-01-01"))


for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_final.RDS")))
  
  #top ten parameter combinations by district
  top_ten_ens <- as.data.frame(tmp[[4]]) %>% mutate(District = district_name, rank = seq(1,10,by = 1)) 
  names(top_ten_ens)[1] <- "num_ens"
  top_ten_df <- rbind(top_ten_df, top_ten_ens)
  
  # rmse for all parameter combinations by district
  rmse <- tmp[[3]] %>% mutate(District = district_name)
  rmse_df <- rbind(rmse_df, rmse)
  
  # individual district level top ten best-fit simulation
  sim_tmp <- tmp[[2]] %>% ungroup() %>% filter(num_ens %in% tmp[[4]]) %>% 
    mutate(District = district_name) %>%
    mutate(total_inc_missing_full = case_when(!is.na(mean) ~ as.numeric(total_inc),
                                         is.na(mean) ~ mean)) %>%
    group_by(num_ens) %>%
    mutate(sim_scaled_full = rescale(total_inc_missing_full)) %>% ungroup()
  sim_df_top_ten <- rbind(sim_df_top_ten, sim_tmp)
  
  
  # top ten paramter combination for each district
  parms_tmp <- parms[tmp[[4]], ] %>% mutate(District = district_name)
  parms_df <- rbind(parms_df, parms_tmp)

}

sim_df_top_ten <- sim_df_top_ten %>% filter(!is.na(District)) %>% filter(year < 2020)    
sim_df_top_ten$date = as.Date(paste(sim_df_top_ten$year, "-", sim_df_top_ten$month, "-", "01", sep = "")) #solve missing date issue for Jessore
top_ten_df <- top_ten_df %>% filter(!is.na(District))  
parms_df <- parms_df %>% filter(!is.na(District))      
rmse_df <- rmse_df %>% filter(!is.na(District))      

save(sim_df_top_ten, top_ten_df, parms_df, rmse_df, file = here("out", "sim_results", "final_sim_results_compiled.RData"))

#### load previously saved data from above ####

load(file = here("out", "sim_results", "final_sim_results_compiled.RData"))

#### top ten parameter combination for each district ####
parms_df_long <- parms_df %>% gather(parameter, value, -c(District))

#### top overall parameter combination and 95% HDI minimizing rmse across all districts ####
best_fit_overall <- rmse_df %>% group_by(num_ens) %>% summarize(sum_rmse = sum(mean_rmse)) %>% arrange(sum_rmse) %>% mutate(rank = seq(1:10000))


parms_top_overall <- parms[best_fit_overall$num_ens[best_fit_overall$rank == 1], ]
parms_top_500_overall <- parms[best_fit_overall$num_ens[best_fit_overall$rank %in% seq(1:500)], ]


#average R0
daily_clim <- readRDS(here("data", "cleaned", "era5_clim_data_for_model.RDS")) %>%
  group_by(date) %>%
  summarize(SH = mean(as.numeric(SH),na.rm=TRUE), tmp = mean(as.numeric(tmp),na.rm=T))

R0 = calc_R0_AH_T_ens(in.parms = parms_top_500_overall, num.ens = 500, sh1 = daily_clim$SH, mean.temp = daily_clim$tmp)[[2]]
parms_top_500_overall$R0 <- colMeans(R0)

R0_top_overall = mean(calc_R0_AH_T_ens(in.parms = parms_top_overall, num.ens = 1, sh1 = daily_clim$SH, mean.temp = daily_clim$tmp)[[2]])

hdi_range <- as.data.frame(cbind(lower = 1, upper = 1, parameter = NA))
varlist =  c("R0","R0max", "R0diff", "qmin", "qmax", "qmid", "Tc", "Texp", "D", "L", "S0perc", "I0")
for (v in 1:length(varlist)){
  lower <- hdi(parms_top_500_overall[,eval(varlist[v])])[1] %>% unlist() %>% as.numeric() %>% round(.,2)
  upper <- hdi(parms_top_500_overall[,eval(varlist[v])])[2] %>% unlist() %>% as.numeric() %>% round(.,2)
  hdi_range <- rbind.data.frame(hdi_range, cbind.data.frame(lower = as.numeric(lower), upper = as.numeric(upper), parameter = varlist[v]))
}

hdi_range <- hdi_range %>% filter(!is.na(parameter))
hdi_range$estimate <- unlist(c(R0_top_overall, parms_top_overall[1:11]))
xtable(hdi_range[,c("parameter", "lower", "upper")])

#### center of mass versus longitude ####

# center of gravity 
sim_df_top_one <- sim_df_top_ten %>% left_join(top_ten_df, by = c("District", "num_ens")) %>% filter(rank == 1)
years <- unique(sim_df_top_one$year)

cm<-cm_sim <- as.data.frame(matrix(NA,ncol=5,nrow=1))
names(cm)<-names(cm_sim) <- c('mean','CIlower','CIupper', 'year', 'District')



for (y in 1:length(years)){
  for(d in 1:length(districts_list_limited)) {
    
    data = sim_df_top_one   %>% mutate(cases = round(mean*100)) %>% filter(District == districts_list_limited[d] & year == years[y]) #observed data
    data_sim = sim_df_top_one   %>% mutate(cases = round(total_inc/1000)) %>% filter(District == districts_list_limited[d] & year == years[y]) #simulated data
    
    if(nrow(data[!is.na(data$cases), ]) == 12 & nrow(data_sim[!is.na(data_sim$cases), ]) == 12) {
      
      cm <- rbind(cm, cbind(center.gravity(data, var = "cases", period = 12), District = districts_list_limited[d], year = years[y]))
      cm_sim <- rbind(cm_sim, cbind(center.gravity(data_sim, var = "cases", period = 12), District = districts_list_limited[d], year = years[y]))
      
    }
  }
  print(y)
}

cm <- cm[-1,]; cm_sim <- cm_sim[-1,]

district_coords <- read.csv(here('data', 'misc', 'sentinel_sites_coords.csv'), header=T) %>%
  group_by(District) %>% summarize(lon = mean(lon, na.rm = TRUE), lat = mean(lat, na.rm = TRUE))

cm_clean <- cm %>% select(year, District, mean) %>%
  left_join(district_coords, by = c("District")) %>% mutate(type = "observed")

cm_sim_clean <- cm_sim %>% select(year, District, mean) %>%
  left_join(district_coords, by = c("District")) %>% mutate(type = "simulated")

cm_combined_allyears <- rbind(cm_clean, cm_sim_clean) 

colors <- c("observed" = "#009E73", "simulated" = "#E69F00")


cm_lon_plt <- cm_combined_allyears %>% #filter(year > 2012) %>%
  ggplot(data = ., aes(x = lon, y = (mean), color = type, shape = type, fill = type)) +
  geom_point(size = 2)+
  #scale_color_manual(values = c('#56B4E9','#E69F00')) +
  #scale_fill_manual(values = c('#56B4E9','#E69F00')) +
  #scale_fill_brewer(palette="Dark2")+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(fill=guide_legend(ncol=2), color = guide_legend(ncol=2)) + 
  geom_smooth(method = "lm") +
  labs(y = "Center of gravity (month)", x = "Longitude ", fill = "", shape = "", color = "") +
  theme_minimal() + 
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))

cm_lat_plt <- cm_combined_allyears %>% #filter(year > 2012) %>%
  ggplot(data = ., aes(x = lat, y = (mean), color = type, shape = type, fill = type)) +
  geom_point(size = 2)+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  #scale_fill_brewer(palette="Dark2")+
  geom_smooth(method = "lm") +
  labs(y = "Center of gravity (month)", x = "Latitude ") +
  theme_minimal() + 
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        legend.position = "none")


# plot seasonality 
sim_seasonality <- sim_df_top_one %>% filter(District %in% districts_list_limited) %>% group_by(month) %>%
  dplyr::summarize(mean = mean(sim_scaled_full, na.rm = TRUE), sd = sd(sim_scaled_full, na.rm = TRUE)) %>% 
  arrange(month) %>% mutate(type = "simulated")

obs_seasonality <- sim_df_top_one %>% filter(District %in% districts_list_limited) %>% 
  group_by(month) %>%
  dplyr::summarize(mean = mean(obs_monthly_scaled_full, na.rm = TRUE), sd = sd(obs_monthly_scaled_full, na.rm = TRUE)) %>% 
  arrange(month) %>% mutate(type = "observed")

seasonality_combined <- rbind(sim_seasonality, obs_seasonality)


seasonality_plt <- seasonality_combined %>% arrange(month) %>%
  mutate(month = as.factor(month)) %>%
  ggplot(aes(x=month,y=mean, color = type, group = type, shape = type)) + 
  geom_point() + theme_minimal() + 
  stat_smooth(method="loess", se=TRUE, aes(fill=type), alpha=0.3) +
  labs(x = "Month") +   coord_cartesian(ylim=c(0,0.53)) +
  scale_y_continuous("Monthly scaled incidence ", breaks = seq(0,0.5,0.1)) + 
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))


#### Fig 4: seasonality/CoM  ####
legend <- get_legend(cm_lon_plt)
plots4 <- plot_grid(seasonality_plt, cm_lat_plt,
                    cm_lon_plt + theme(legend.position = "none"), 
                    nrow=1, rel_widths = c(1,1,1),
                    labels = c("A","B", "C"))

fig_sim_results <-  plot_grid(plot_grid(plots4, plot_grid(legend,NULL, ncol=1), 
                             ncol = 1,rel_heights=c(1, 0.1)), ncol = 1,
                             NULL, rel_heights=c(1, 0.02))                                               
fig_sim_results
ggsave(fig_sim_results, file = here("out", "plots", "Fig4_fig_sim_results_v3.png"), width = 7.5, height = 3)

#### Fig 3: best-fit simulation plots (top ten by district) ####

plot_districts <- c("Barisal", "Sylhet", "Chittagong","Dhaka")
colors <- c("observed" = "#009E73", "simulated" = "#E69F00")
linetype = c("observed" = 2, "simulated" = 1)

sim_plot <- sim_df_top_ten %>% mutate(date = as_date(date)) %>%
  filter(District %in% plot_districts) %>%
  ggplot() + 
  annotate(geom= 'rect', xmin=as.Date("2018-01-01"),xmax=as.Date("2019-12-31"),
           ymin=0,ymax=1, fill = "grey", alpha=0.3)+
  geom_line(data = sim_df_top_one %>% filter(District %in% plot_districts), 
             aes(x = date, y = obs_monthly_scaled_full, color = 'observed', lty = 'observed'), size = 1) +
  geom_line(aes(x = date , y = sim_scaled_full, group = as.factor(num_ens), color = 'simulated', lty = 'simulated'), alpha = 0.3) +
  theme_minimal() + guides(linetype = "none") + 
  scale_color_manual(values = colors)+ 
  scale_linetype_manual(values = linetype) + 
  labs(color = "") +   facet_wrap(~District, scales = "free") + 
  scale_x_date("",limits = c(as.Date("2012-01-01"), as.Date("2019-12-31"))) + 
  scale_y_continuous("Monthly scaled incidence", limits = c(0,1)) + 
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
sim_plot
ggsave(sim_plot, file  = here("out", "plots", "Fig3_simulations_4_districts_v2.png"), 
       width = 7, height = 5)



sim_df_top_ten_max_min <- sim_df_top_ten %>% group_by(date, District, obs_monthly_scaled_full) %>% 
  filter(!is.na(sim_scaled_full)) %>%
  summarize(sim_scaled_max = max(sim_scaled_full, na.rm = T), 
            sim_scaled_min = min(sim_scaled_full, na.rm = T)) %>%
  mutate(in_bound = case_when(round(obs_monthly_scaled_full,1) >= round(sim_scaled_min,1) & round(obs_monthly_scaled_full,1) <= round(sim_scaled_max,1) ~ 1,
                              TRUE ~ 0))
sum(sim_df_top_ten_max_min$in_bound) / nrow(sim_df_top_ten_max_min) #% in bound

sim_plot <- sim_df_top_ten %>% mutate(date = as_date(date)) %>%
  filter(District %in% plot_districts) %>%
  ggplot() + 
  annotate(geom= 'rect', xmin=as.Date("2018-01-01"),xmax=as.Date("2019-12-31"),
           ymin=0,ymax=1, fill = "grey", alpha=0.3)+
  geom_line(data = sim_df_top_one %>% filter(District %in% plot_districts), 
            aes(x = date, y = obs_monthly_scaled_full, color = 'observed', lty = 'observed'), size = 1) +
  geom_line(data = sim_df_top_one %>% filter(District %in% plot_districts), 
            aes(x = date, y = sim_scaled_full, color = 'simulated', lty = 'simulated'), size = 0.8) +
  #
  #geom_line(data = sim_df_top_ten_max_min %>% filter(District %in% plot_districts), 
  #          aes(x = date, y = sim_scaled_max, color = 'simulated', lty = 'simulated'), alpha = 0.5) +
  #geom_line(data = sim_df_top_ten_max_min %>% filter(District %in% plot_districts), 
  #          aes(x = date, y = sim_scaled_min, color = 'simulated', lty = 'simulated'), alpha = 0.5) +
  
  geom_ribbon(data = sim_df_top_ten_max_min %>% filter(District %in% plot_districts), 
              aes(x = date,
                  ymin = sim_scaled_min,
                  ymax = sim_scaled_max), alpha = 0.5, fill = "#E69F00") +
  
  #geom_line(aes(x = date , y = sim_scaled, group = as.factor(num_ens), color = 'simulated', lty = 'simulated'), alpha = 0.3) +
  theme_minimal() + guides(linetype = "none") + 
  scale_color_manual(values = colors)+ 
  scale_linetype_manual(values = linetype) + 
  labs(color = "") +   facet_wrap(~District, scales = "free") + 
  scale_x_date("",limits = c(as.Date("2012-01-01"), as.Date("2019-12-31"))) + 
  scale_y_continuous("Monthly scaled incidence", limits = c(0,1)) + 
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
sim_plot
ggsave(sim_plot, file  = here("out", "plots", "Fig3_simulations_4_districts_v3.png"), 
       width = 7, height = 5)


# plot all districts for supplementary material
sim_plot_all_districts <- sim_df_top_ten %>% mutate(date = as_date(date)) %>%
  #filter(District %in% plot_districts) %>%
  ggplot() + 
  annotate(geom= 'rect', xmin=as.Date("2018-01-01"),xmax=as.Date("2019-12-31"),
           ymin=0,ymax=1, fill = "grey", alpha=0.3)+
  geom_line(data = sim_df_top_one, #%>% filter(District %in% plot_districts), 
            aes(x = date, y = obs_monthly_scaled_full, color = 'observed', lty = 'observed'), size = 1) +
  geom_line(aes(x = date , y = sim_scaled_full, group = as.factor(num_ens), color = 'simulated', lty = 'simulated'), alpha = 0.3) +
  theme_minimal() + guides(linetype = "none") + 
  scale_color_manual(values = colors)+ 
  scale_linetype_manual(values = linetype) + 
  labs(color = "") +   facet_wrap(~District, scales = "free", ncol = 3) + 
  scale_x_date("",limits = c(as.Date("2012-01-01"), as.Date("2019-12-31"))) + 
  scale_y_continuous("Monthly scaled incidence", limits = c(0,1)) + 
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
sim_plot_all_districts
ggsave(sim_plot_all_districts, file  = here("out", "plots", "simulations_all_districts.png"), width = 6, height = 6)



#### climate change projections ####

plot_districts <- c("Barisal", "Sylhet", "Chittagong","Dhaka")
proj_source = "CMIP6_GFDL"

# get best-fit simulation and climate change simulation
num_sim = 1
sim_df_cc <- as.data.frame(cbind(inc = NA, District = NA, type = NA, month = NA, year = NA)) %>% mutate(date = as.Date("2012-01-01"))
for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_final.RDS")))
  sim_tmp <- tmp[[5]] %>% ungroup() %>% filter(num_ens == tmp[[4]][1]) %>% 
    select(-num_ens) %>% mutate(District = district_name, type = "baseline")
  sim_tmp_cc <-   readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_projection_",proj_source,".RDS")))[[2]]  %>% 
    mutate(District = district_name, type = "projection")

  
  sim_df_cc <- rbind(sim_df_cc, sim_tmp, sim_tmp_cc)
  rm(sim_tmp, sim_tmp_cc, tmp)
  
  
}


sim_df_cc <- sim_df_cc %>% filter(!is.na(District))

clim <- readRDS(file = here("data", "cleaned", paste0("era5_",proj_source,"_modified_clim_data_for_projection.RDS"))) %>%
  mutate(week = lubridate::week(date)) 


#### Fig 5: climate change plot ####

cc_seasonality_plot <- sim_df_cc %>% 
  filter(District %in% plot_districts) %>%
  mutate(doy = strftime(date, format = "%j"), 
                     week = lubridate::week(date)) %>% group_by(week, District, type) %>%
  dplyr::summarize(mean = mean(inc, na.rm = TRUE), 
                   sd = sd(inc, na.rm = TRUE)) %>% arrange(week) %>%
  ggplot(.,aes(x = week , y = mean, color = type)) + geom_line(size=1) +
  #scale_color_manual(values=wes_palette(n=2, name="Darjeeling1"))+
  scale_colour_manual("", values = c("coral", "cyan4")) +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  facet_wrap(~District, scales = "free", ncol = 2) +
  theme_minimal() + labs(y = "Mean incidence") +
  scale_x_continuous("Week", breaks = seq(0,52,10)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12),
        legend.position = "bottom")

legend <- get_legend(cc_seasonality_plot)
#ggsave(cc_seasonality_plot, file  = here("out", "plots", "cc_seasonality_4_districts.png"), width = 6, height = 5)

SH_cc_plot <- clim %>% group_by(week) %>%
  summarize(SH  = mean(SH, na.rm = TRUE), SH_cc = mean(SH_cc, na.rm = TRUE)) %>%
  ggplot() + theme_minimal() +
  geom_point(aes(x = week, y = SH), color = "coral") +
  geom_point(aes(x = week, y = SH_cc), color = "cyan4") +
  labs(y = "SH (g/kg)") +
  scale_x_continuous("Week", breaks = seq(0,52,10)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))

tmp_cc_plot <- clim %>% group_by(week) %>%
  summarize(tmp  = mean(tmp, na.rm = TRUE), tmp_cc = mean(tmp_cc, na.rm = TRUE)) %>%
  ggplot() + geom_point(aes(x = week, y = tmp), color = "coral") +
  geom_point(aes(x = week, y = tmp_cc), color = "cyan4") +
  labs(y = "Temperature (C)    ") + theme_minimal() +
  scale_x_continuous("Week", breaks = seq(0,52,10)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))



cc_combined <- ggarrange(ggarrange(
  plot_grid(SH_cc_plot, tmp_cc_plot, ncol = 1, labels = c("A", "B")),
  cc_seasonality_plot + theme(legend.position = "none"),
  ncol = 2, labels = c("", "C"), widths = c(0.5,1)),
  legend, nrow = 2, heights = c(1,0.10))


cc_combined
ggsave(cc_combined, width = 9, height = 4.5,
       file = here("out", "plots", paste0("Fig5_cc_combined_",proj_source,"_v2.png")))

# plot cc seasonality for all districts
cc_seasonality_plot_all_districts <- sim_df_cc %>% 
  mutate(doy = strftime(date, format = "%j"), 
         week = lubridate::week(date)) %>% group_by(week, District, type) %>%
  dplyr::summarize(mean = mean(inc, na.rm = TRUE), 
                   sd = sd(inc, na.rm = TRUE)) %>% arrange(week) %>%
  ggplot(.,aes(x = week , y = mean, color = type)) + geom_line(size=1) +
  #scale_color_manual(values=wes_palette(n=2, name="Darjeeling1"))+
  scale_colour_manual("", values = c("coral", "cyan4")) +
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  facet_wrap(~District, scales = "free") +
  theme_minimal() + labs(y = "Mean incidence") +
  scale_x_continuous("Week", breaks = seq(0,52,10)) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12),
        legend.position = "bottom")
cc_seasonality_plot_all_districts

ggsave(cc_seasonality_plot_all_districts, file  = here("out", "plots", "cc_seasonality_plot_all_districts.png"), width = 6, height = 5)


#### SM Fig: RMSE  ####
rmse_df <- sim_df_top_ten %>% mutate(period = case_when(year <= 2017 ~ "Training",
                                                                  TRUE ~ "Testing")) %>% 
  mutate(error = (obs_monthly_scaled_full - sim_scaled_full)^2) %>%
  group_by(District, num_ens, period) %>% 
  summarize(sse = sum(error, na.rm = TRUE), n = length(!is.na(obs_monthly_scaled_full)),  .groups = 'drop') %>%
  mutate(rmse = sqrt(sse/n)) 

rmse_df$rmse[rmse_df$rmse == 0 & rmse_df$period == "Testing" & rmse_df$District %in% c("Bogra" , "Mymensingh")] <- NA


rmse_plot <- ggplot(rmse_df , aes(x = factor(period), y = rmse)) +
  geom_boxplot() +
  labs(color = "", x = "", y = "RMSE")+
  scale_y_continuous("Root Mean Squared Error", limits = c(0, 0.42)) + 
  scale_x_discrete("") + 
  theme_minimal() +  facet_wrap(~District, ncol = 4, scales = "free") +
  theme(legend.position = "bottom", 
                        axis.text = element_text(color = "black"),
                        axis.ticks = element_line(color = "black"),
                        axis.line = element_line(color = "black", size = 0.3),
                        strip.text = element_text(color = "black", size = 12))
rmse_plot
ggsave(rmse_plot, file  = here("out", "plots", "rmse_all_districts.png"), width = 7, height = 7)
rmse_df%>% summarize(mean(rmse, na.rm = TRUE))


#### RMSE for NULL 1 ####
num_sim = 1
sim_df_top_one_NULL1<- as.data.frame(cbind(sim_scaled = NA, num_ens = NA, year = NA, month = NA, District = NA, total_inc = NA, mean = NA, mean_test_data = NA, 
                                           obs_monthly_scaled = NA, obs_monthly_scaled_full = NA, total_inc_missing = NA,
                                           total_inc_missing_full = NA, sim_scaled_full = NA)) %>% mutate(date = as.Date("2012-01-01"))


for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_final_NULL1.RDS")))
  
  # individual district level top ten best-fit simulation
  sim_tmp <- tmp[[2]] %>% ungroup() %>% filter(num_ens %in% tmp[[4]][1]) %>% 
    mutate(District = district_name) %>%
    mutate(total_inc_missing_full = case_when(!is.na(mean) ~ as.numeric(total_inc),
                                              is.na(mean) ~ mean)) %>%
    group_by(num_ens) %>%
    mutate(sim_scaled_full = rescale(total_inc_missing_full)) %>% ungroup()
  sim_df_top_one_NULL1 <- rbind(sim_df_top_one_NULL1, sim_tmp)
  
 
}

sim_df_top_one_NULL1 <- sim_df_top_one_NULL1 %>% filter(!is.na(District)) %>% filter(year < 2020)    
sim_df_top_one_NULL1$date = as.Date(paste(sim_df_top_one_NULL1$year, "-", sim_df_top_one_NULL1$month, "-", "01", sep = "")) #solve missing date issue for Jessore

rmse_df_NULL1 <- sim_df_top_one_NULL1 %>% mutate(period = case_when(year <= 2017 ~ "Training",
                                                        TRUE ~ "Testing")) %>% 
  mutate(error = (obs_monthly_scaled_full - sim_scaled_full)^2) %>%
  group_by(District, num_ens, period) %>% 
  summarize(sse = sum(error, na.rm = TRUE), n = length(!is.na(obs_monthly_scaled_full)),  .groups = 'drop') %>%
  mutate(rmse = sqrt(sse/n)) 

rmse_df_NULL1$rmse[rmse_df_NULL1$rmse == 0 & rmse_df_NULL1$period == "Testing" & rmse_df_NULL1$District %in% c("Bogra" , "Mymensingh")] <- NA
rmse_df_NULL1 %>% summarize(mean(rmse, na.rm = TRUE))


#### RMSE for NULL 2 ####
num_sim = 1
rmse_all_NULL2 <- as.data.frame(cbind(mean_rmse= NA, num_ens = NA, District = NA))
sim_df_top_one_NULL2<- as.data.frame(cbind(sim_scaled = NA, num_ens = NA, year = NA, month = NA, District = NA, total_inc = NA, mean = NA, mean_test_data = NA, 
                                           obs_monthly_scaled = NA, obs_monthly_scaled_full = NA, total_inc_missing = NA,
                                           total_inc_missing_full = NA, sim_scaled_full = NA)) %>% mutate(date = as.Date("2012-01-01"))

for (d in 1:length(districts_list_limited)) {
  
  # set name of district
  district_name = districts_list_limited[d]
  tmp <- readRDS(file = here("out", "sim_results", paste0(district_name, "_numsim_",num_sim,"_topten_final_NULL2.RDS")))
  
  
  # rmse for all parameter combinations by district
  rmse_tmp <- tmp[[3]] %>% mutate(District = district_name)
  rmse_all_NULL2 <- rbind(rmse_all_NULL2, rmse_tmp)
  
  
  # individual district level top ten best-fit simulation
  sim_tmp <- tmp[[2]] %>% ungroup() %>% filter(num_ens %in% tmp[[4]][1]) %>% 
    mutate(District = district_name) %>%
    mutate(total_inc_missing_full = case_when(!is.na(mean) ~ as.numeric(total_inc),
                                              is.na(mean) ~ mean)) %>%
    group_by(num_ens) %>%
    mutate(sim_scaled_full = rescale(total_inc_missing_full)) %>% ungroup()
  sim_df_top_one_NULL2 <- rbind(sim_df_top_one_NULL2, sim_tmp)
  
  
}

sim_df_top_one_NULL2 <- sim_df_top_one_NULL2 %>% filter(!is.na(District)) %>% filter(year < 2020)  
sim_df_top_one_NULL2$date = as.Date(paste(sim_df_top_one_NULL2$year, "-", sim_df_top_one_NULL2$month, "-", "01", sep = "")) #solve missing date issue for Jessore

rmse_all_NULL2 <- rmse_all_NULL2 %>% filter(!is.na(District))
# rmse for best-fit simulation
rmse_df_NULL2 <- sim_df_top_one_NULL2 %>% mutate(period = case_when(year <= 2017 ~ "Training",
                                                                    TRUE ~ "Testing")) %>% 
  mutate(error = (obs_monthly_scaled_full - sim_scaled_full)^2) %>%
  group_by(District, num_ens, period) %>% 
  summarize(sse = sum(error, na.rm = TRUE), n = length(!is.na(obs_monthly_scaled_full)),  .groups = 'drop') %>%
  mutate(rmse = sqrt(sse/n)) 

rmse_df_NULL2$rmse[rmse_df_NULL2$rmse == 0 & rmse_df_NULL2$period == "Testing" & rmse_df_NULL2$District %in% c("Bogra" , "Mymensingh")] <- NA
rmse_df_NULL2 %>% summarize(mean(rmse, na.rm = TRUE))


#### SM Fig: RMSE all ####
rmse_df$type = "Climate forcing"
rmse_df_NULL1$type = "No forcing"
rmse_df_NULL2$type = "Sinusoidal forcing"

rmse_df_top_one <- rmse_df %>% left_join(top_ten_df, by = c("District", "num_ens")) %>% 
  filter(rank == 1) %>% select(-rank)

rmse_df_all <- rbind(rmse_df_NULL1 , 
                     rmse_df_NULL2 , 
                     rmse_df_top_one ) %>%
  mutate(rmse_clean = case_when (is.infinite(rmse) ~ as.numeric(NA),
                                 TRUE ~ as.numeric(rmse)))

rmse_df_all$type <- factor(rmse_df_all$type, levels=c('No forcing', 'Sinusoidal forcing', 'Climate forcing'))

rmse_plot_all <- ggplot(rmse_df_all , aes(x = factor(period), y = rmse_clean, col = type)) +
  geom_point(position = position_dodge(width = 0.3), cex = 2.5) +
  labs(color = "", x = "", y = "RMSE")+
  scale_y_continuous("Root Mean Squared Error", limits = c(0,0.65)) + 
  scale_x_discrete("") + 
  theme_minimal() +  facet_wrap(~District, ncol = 4, scales = "free") +
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
rmse_plot_all
ggsave(rmse_plot_all, file  = here("out", "plots", "rmse_all_districts_all.png"), width = 7, height = 7)


rmse_df_all %>% group_by(type) %>% summarize(mean_rmse = mean(rmse, na.rm = TRUE))


#### NULL 2 R0 ####
parms_NULL2 <- readRDS(here("data", "cleaned", "LHS_10000_parameter_matrix_tuning_round4_NULL2.RDS"))

best_fit_overall_NULL2 <- rmse_all_NULL2 %>% group_by(num_ens) %>% summarize(sum_rmse = sum(mean_rmse)) %>% arrange(sum_rmse) %>% mutate(rank = seq(1:10000))
parms_top_overall_NULL2 <- parms_NULL2[best_fit_overall_NULL2$num_ens[best_fit_overall_NULL2$rank == 1], ]

R0_NULL2 <- calc_R0_sinusoidal_ens(in.parms = parms_top_overall_NULL2, num.ens = 1, sh1 = seq(1:365))
R0_NULL2_df <- as.data.frame(cbind(R0 = R0_NULL2, doy = seq(1:365)))

names(R0_NULL2_df)[1] <- "R0"

R0_NULL2_plot <- ggplot(R0_NULL2_df, aes(x=doy, y=R0)) +
  #geom_point() +
  geom_line() +
  theme_minimal()+labs(x = "Day of year", y = parse(text = "R[0]"))+
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))

ggsave(R0_NULL2_plot, file  = here("out", "plots", "R0_NULL2_plot.png"), width = 4, height = 4)

#### SM fig: compare with NULL2 ####
colors <- c("observed" = "#009E73", "climate forcing" = "#E69F00", "sinusoidal forcing" = "light blue")
linetype = c("observed" = 2, "climate forcing" = 1, "sinusoidal forcing" = 1)

sim_plot_compare_NULL2 <- sim_df_top_ten %>% mutate(date = as_date(date)) %>%
  ggplot() + 
  annotate(geom= 'rect', xmin=as.Date("2018-01-01"),xmax=as.Date("2019-12-31"),
           ymin=0,ymax=1, fill = "grey", alpha=0.3)+
  geom_line(data = sim_df_top_one_NULL2 %>% filter(!is.na(District)), #%>% filter(District %in% plot_districts), 
            aes(x = date, y = sim_scaled_full, color = 'sinusoidal forcing', lty = 'sinusoidal forcing'), size = 1) +
  
  geom_line(data = sim_df_top_one, #%>% filter(District %in% plot_districts), 
            aes(x = date, y = obs_monthly_scaled_full, color = 'observed', lty = 'observed'), size = 1) +
  geom_line(aes(x = date , y = sim_scaled_full, group = as.factor(num_ens), color = 'climate forcing', lty = 'climate forcing'), alpha = 0.3) +
  theme_minimal() + guides(linetype = "none") + 
  scale_color_manual(values = colors)+ 
  scale_linetype_manual(values = linetype) + 
  labs(color = "") +   facet_wrap(~District, scales = "free", ncol = 3) + 
  scale_x_date("",limits = c(as.Date("2012-01-01"), as.Date("2019-12-31"))) + 
  scale_y_continuous("Monthly scaled incidence", limits = c(0,1)) + 
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
sim_plot_compare_NULL2
ggsave(sim_plot_compare_NULL2, file  = here("out", "plots", "sim_plot_compare_NULL2.png"), width = 7, height = 6)


# split into two plots

sim_plot_compare_NULL2_1 <- sim_df_top_ten %>% mutate(date = as_date(date)) %>%
  filter(District %in% districts_list_limited[1:6]) %>%
  ggplot() + 
  annotate(geom= 'rect', xmin=as.Date("2018-01-01"),xmax=as.Date("2019-12-31"),
           ymin=0,ymax=1, fill = "grey", alpha=0.3)+
  geom_line(data = sim_df_top_one_NULL2 %>% filter(District %in% districts_list_limited[1:6]), #%>% filter(District %in% plot_districts), 
            aes(x = date, y = sim_scaled_full, color = 'sinusoidal forcing', lty = 'sinusoidal forcing'), size = 1, na.rm = TRUE) +
  
  geom_line(data = sim_df_top_one %>% filter(District %in% districts_list_limited[1:6]), 
            aes(x = date, y = obs_monthly_scaled_full, color = 'observed', lty = 'observed'), size = 1, na.rm = TRUE) +
  geom_line(aes(x = date , y = sim_scaled_full, group = as.factor(num_ens), color = 'climate forcing', lty = 'climate forcing'), alpha = 0.3, na.rm = TRUE) +
  theme_minimal() + guides(linetype = "none") + 
  scale_color_manual(values = colors)+ 
  scale_linetype_manual(values = linetype) + 
  labs(color = "") +   facet_wrap(~District, scales = "free", ncol = 2) + 
  scale_x_date("",limits = c(as.Date("2012-01-01"), as.Date("2019-12-31"))) + 
  scale_y_continuous("Monthly scaled incidence", limits = c(0,1)) + 
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
sim_plot_compare_NULL2_1
ggsave(sim_plot_compare_NULL2_1, file  = here("out", "plots", "sim_plot_compare_NULL2_1.png"), width = 7, height = 6)


sim_plot_compare_NULL2_2 <- sim_df_top_ten %>% mutate(date = as_date(date)) %>%
  filter(District %in% districts_list_limited[7:12]) %>%
  ggplot() + 
  annotate(geom= 'rect', xmin=as.Date("2018-01-01"),xmax=as.Date("2019-12-31"),
           ymin=0,ymax=1, fill = "grey", alpha=0.3)+
  geom_line(data = sim_df_top_one_NULL2 %>% filter(District %in% districts_list_limited[7:12]), #%>% filter(District %in% plot_districts), 
            aes(x = date, y = sim_scaled_full, color = 'sinusoidal forcing', lty = 'sinusoidal forcing'), size = 1, na.rm = TRUE) +
  
  geom_line(data = sim_df_top_one %>% filter(District %in% districts_list_limited[7:12]), 
            aes(x = date, y = obs_monthly_scaled_full, color = 'observed', lty = 'observed'), size = 1, na.rm = TRUE) +
  geom_line(aes(x = date , y = sim_scaled_full, group = as.factor(num_ens), color = 'climate forcing', lty = 'climate forcing'), alpha = 0.3, na.rm = TRUE) +
  theme_minimal() + guides(linetype = "none") + 
  scale_color_manual(values = colors)+ 
  scale_linetype_manual(values = linetype) + 
  labs(color = "") +   facet_wrap(~District, scales = "free", ncol = 2) + 
  scale_x_date("",limits = c(as.Date("2012-01-01"), as.Date("2019-12-31"))) + 
  scale_y_continuous("Monthly scaled incidence", limits = c(0,1)) + 
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
sim_plot_compare_NULL2_2
ggsave(sim_plot_compare_NULL2_2, file  = here("out", "plots", "sim_plot_compare_NULL2_2.png"), width = 7, height = 6)



#### SM fig: variation across districts for top ten parameter estimates ####
label <- as_labeller(c(D = "D", I0 = "I[0]", L = "L", qmax = "q[max]",
                       qmin = "q[min]", qmid = "q[mid]", R0diff = "R['0,diff']", 
                       R0max = "R['0,max']", S0perc = "S['0']~(percent)", Tc = "T[c]", 
                       Texp = "T[exp]"),
                     default = label_parsed)

parms_top_ten_district <- ggplot(parms_df_long %>% filter(parameter != "Tdiff"), aes(x = parameter, y = value)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_x_discrete("")+
  theme_minimal() +
  facet_wrap(~parameter, labeller = label, scales = "free", ncol = 3) +
  theme(legend.position = "bottom", 
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.line.y = element_line(color = "black", size = 0.3),
        strip.text.y = element_text(color = "black", size = 12),
        axis.text.x=element_blank())
parms_top_ten_district
ggsave(parms_top_ten_district, file  = here("out", "plots", "parms_top_ten_district.png"), width = 6, height = 6)

#### SM fig: relationship between R0 and climate ####
R0 = calc_R0_AH_T_ens(in.parms = parms_top_overall, num.ens = 1, sh1 = seq(0,25), mean.temp = rep(20,26))[[2]]
R01 = calc_R0_AH_T_ens(in.parms = parms_top_overall, num.ens = 1, sh1 = seq(0,25), mean.temp = rep(25,26))[[2]]
R02 = calc_R0_AH_T_ens(in.parms = parms_top_overall, num.ens = 1, sh1 = seq(0,25), mean.temp = rep(27,26))[[2]]
R03 = calc_R0_AH_T_ens(in.parms = parms_top_overall, num.ens = 1, sh1 = seq(0,25), mean.temp = rep(30,26))[[2]]

R0df_plot <- rbind(cbind(R = R0, SH = seq(0,25), temp = 20), cbind(R = R01, SH = seq(0,25), temp = 25),
                   cbind(R = R02, SH = seq(0,25), temp = 27), cbind(R=R03, SH = seq(0,25), temp = 29)) %>% as.data.frame()
names(R0df_plot)[1] <- "R0"

R0_clim_plot <- ggplot(R0df_plot, aes(x=SH, y=R0, group = factor(temp), col = factor(temp))) +
  #geom_point() +
  geom_line() +
  scale_color_brewer("Temperature",palette = "Blues")+
  theme_minimal()+labs(y = parse(text = "R[0]"))+
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12))
  
ggsave(R0_clim_plot, file  = here("out", "plots", "R0_clim_plot.png"), width = 4, height = 4)

