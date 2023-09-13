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
library(texreg)
library(car)
library(mgcv)

rm(list = ls())
source(here("code", "cg_circular.R"))

## Notes:
# Assuming NISB and HBIS data are extracted and cleaned and saved in data > cleaned folder
# Assuming precipitation from CRU data are extracted and cleaned and saved in data > cleaned > clim_all.rds 



#### epi data (from two sources) ####

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


#### climate data ####
tmp_dat <- readRDS(here('data', 'cleaned', 'all_tmp_by_site_era5.rds')) %>%
  distinct(month, year, Hospital, District, monthlytmp, monthlytmpK) %>%
  rename(tmp = monthlytmp, tmp_K = monthlytmpK)

pre_dat <- readRDS(file = here('data', 'cleaned', 'clim_all.rds')) %>%
  distinct(month, year, Hospital, District, pre)

clim_dat <- readRDS(here('data', 'cleaned', 'all_SH_by_site_era5.rds')) %>%
  distinct(month, year, Hospital, District, monthlySH) %>%
  left_join(.,tmp_dat, by = c("month", "year", "Hospital", "District")) %>%
  left_join(., pre_dat, by = c("month", "year", "Hospital", "District")) %>%
  rename(SH = monthlySH)


#### merge case data and climate data ####

df_flu_clim <- df_flu_all %>% left_join(., clim_dat, by = c("month", "year", "Hospital", "District")) %>%
  filter(!is.na(tmp)) %>%
  mutate(season = case_when(month %in% c(1,2,12) ~ "Winter",
                            month %in% c(3:5) ~ "Pre-Monsoon",
                            month %in% c(6:9) ~ "Monsoon",
                            TRUE ~ "Post-Monsoon")) %>%
  mutate(Positive = as.numeric(Flu),
         Sampled = as.numeric(Sample),
         Prop = as.numeric(Prop)) %>%
  filter(Sampled > 0) %>% 
  mutate(logSampled = log(Sampled), season = as.factor(season), year = as.factor(year), Hospital = as.factor(Hospital)) 




#### plot number sampled and number tested positive (supplementary figure 1) ####
colors <- c("Cases" = "palegreen3", "Total sampled" = "red", "Proportion positive" = "Dodgerblue4")

cases_sampled_plot <- df_flu_all %>% filter(as.numeric(year) < 2020) %>% group_by(date, month, year) %>% summarize(Sample = sum(Sample), Flu = sum(Flu) ) %>%
  ggplot(.) + geom_line(aes(x = date, y = Flu, col = "Cases"), cex = 1.5) +
  geom_line(aes(x = date, y = Sample, col = "Total sampled"), linetype = "dashed") +
  geom_line(aes(x = date, y = 1000*Flu/Sample, col = "Proportion positive")) +
  scale_color_manual(values = colors)+ 
  scale_fill_manual(values = colors)+ 
  scale_y_continuous(sec.axis=sec_axis(~.*0.001,name="Proportion positive for influenza")) +
  labs(y = "Total sampled and Influenza-positive cases", x = " ", color = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key = element_rect(colour = NA, fill = NA),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))
cases_sampled_plot

ggsave(cases_sampled_plot, file = here("out", "plots", "cases_sampled_plot.png"), width = 8, height = 7)


#### regression analysis ####

df_flu_clim <- df_flu_clim %>% as_tibble %>% mutate(month = as.factor(month))
fit1 <- glm(Prop ~ tmp   + I(tmp^2) + SH  + I(SH^2) + pre + I(pre^2) + month + 
              Hospital, data = df_flu_clim %>% filter(as.numeric(year) < 2020), 
            family = binomial, weights = Sampled)
summary(fit1)

fit2 <- glm(Prop ~ tmp   + I(tmp^2)  + as.factor(month) + Hospital, data = df_flu_clim %>% filter(as.numeric(year) < 2020), family = binomial, weights = Sampled)
fit3 <- glm(Prop ~ SH  + I(SH^2) +  as.factor(month) + Hospital, data = df_flu_clim %>% filter(as.numeric(year) < 2020), family = binomial, weights = Sampled)
fit4 <- glm(Prop ~ pre + I(pre^2) + as.factor(month) + Hospital, data = df_flu_clim %>% filter(as.numeric(year) < 2020), family = binomial, weights = Sampled)

# jointly test significance of SH and SH^2
linearHypothesis(fit1, c("SH=0", "I(SH^2)=0"), white.adjust = "hc1")



model_results1 <- get_model_data(fit1, type = "est", robust = TRUE, show.intercept = TRUE)
model_results2 <- get_model_data(fit2, type = "est", robust = TRUE, show.intercept = TRUE)
model_results3 <- get_model_data(fit3, type = "est", robust = TRUE, show.intercept = TRUE)
model_results4 <- get_model_data(fit4, type = "est", robust = TRUE, show.intercept = TRUE)


# regression table
tex1 <- texreg::extract(fit1)
tex1@coef <- model_results1$estimate
tex1@ci.low <- model_results1$conf.low
tex1@ci.up <- model_results1$conf.high


tex2 <- texreg::extract(fit2)
tex2@coef <- model_results2$estimate
tex2@ci.low <- model_results2$conf.low
tex2@ci.up <- model_results2$conf.high


tex3 <- texreg::extract(fit3)
tex3@coef <- model_results3$estimate
tex3@ci.low <- model_results3$conf.low
tex3@ci.up <- model_results3$conf.high


tex4 <- texreg::extract(fit4)
tex4@coef <- model_results4$estimate
tex4@ci.low <- model_results4$conf.low
tex4@ci.up <- model_results4$conf.high

texreg(c(tex1,tex2, tex3, tex4), ci.force = TRUE, ci.test = 1,
       doctype=FALSE,  align.center=TRUE,
       caption= "")


##### Figure 1 #####

climate_season_plot <- df_flu_clim %>% arrange(Prop) %>% 
  mutate(season = fct_relevel(season, "Pre-Monsoon")) %>%
  ggplot() + facet_wrap(~season, nrow = 1, scales = "free_y") + 
  geom_point(aes(x=tmp, y = SH, color = Prop), alpha = 0.5, size = 2) +
  scale_y_continuous("Specific Humidity (g/kg)", limits = c(7,21.4),
                     breaks = seq(from=8,to=22,by=4)) + 
  scale_x_continuous("Temperature (Celsius)", limits = c(17.15,32),
                     breaks = seq(from=18,to=32,by=4)) + 
  scale_color_viridis(option = "B", direction = -1, 
                      "Proportion positive for influenza\n") +
  theme_minimal() + guides(size = "none") + 
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3),
        strip.text = element_text(color = "black", size = 12),
        legend.position = "bottom") 
climate_season_plot
ggsave(climate_season_plot, file = here("out", "plots", "Fig1_climate_season_proportion_v2.png"), width = 8, height = 3)


##### Figure 2 ####
# plot predicted probabilities by covariate
precip_plot <- plot_model(fit1, type = "pred", terms = c("pre[all]"),title = "",
                          axis.title = c("Precipitation (mm)", "Predicted proportion \npositive for influenza")) +
  theme_minimal() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))

SH_plot <- plot_model(fit1, type = "pred", terms = c("SH[0:23]"), 
                      title = "", 
                      axis.title = c("SH(g/kg)", "Predicted proportion \npositive for influenza")) + 
  theme_minimal() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))
SH_plot


tmean_plot <- plot_model(fit1, type = "pred", terms = c("tmp[all]"),
                         title = "", 
                         axis.title = c("Temperature (Celsius)", "Predicted proportion \npositive for influenza")) + 
  theme_minimal() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))
tmean_plot


season_plot <- plot_model(fit1, type = "pred", terms =c("month"), title = "", 
                          axis.title = c("Month", "Predicted proportion \npositive for influenza")) + 
  theme_minimal() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black", size = 0.3))
season_plot

reg_plots <- plot_grid(SH_plot, tmean_plot, precip_plot, season_plot, labels = c("A", "B", "C", "D"))
ggsave(reg_plots, file = here("out", "plots", "Fig2_all_conditional_effects_v2.png"), width = 6, height = 6)


