################################################################################
#                             Reproducibility Code                             #
#       Paper: Cellwise outlier detection in heterogeneous populations         #
#   Authors: Zaccaria, G., García-Escudero, L.A., Greselin, F., Mayo-Íscar, A. #
#          Simulation Study - Main Article - Section 3, Figures 1 and 2        #
################################################################################

## To reproduce an entire figure from the simulation study, as reported 
## in the Main Article or in the Supplementary Material, it is necessary to run 
## all scenarios for each level of contamination or at least one scenario of the 
## simulation study for all three level of contamination (0%, 5%, and 10%).
## We have attach an .RData file containing the results of the three scenarios reported 
## in the Main Article, called "Data for Figure Reproducibility.RData", which can 
## be used to fully reproduce Figures 1 and 2.

## PLEASE, BEFORE RUNNING THE CODE READ CAREFULLY THE INSTRUCTIONS REPORTED IN THE
## README.md (OR .html) FILE AND INSTALL THE PACKAGES IF NECESSARY.

## Set the directory where the folder named "Code for Technometrics Reproducibility" 
## currently is.
setwd("path/Code for Technometrics Reproducibility")

## Load the required packages and the data.
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

load("Data for Figure Reproducibility.RData")

## Prepare the object for ggplot.
df_long <- stat_ggplot %>%
  pivot_longer(cols = c(MR, MSEmu1, MSEmu2, KLsigma1, KLsigma2), 
               names_to = "Metric", 
               values_to = "Value")

df_long_Scen3 <- stat_ggplotscen3 %>%
  pivot_longer(cols = c(MR, MSEmu1, MSEmu2, MSEmu3, MSEmu4, KLsigma1, KLsigma2, KLsigma3, KLsigma4), 
               names_to = "Metric", 
               values_to = "Value")

df_tot <- rbind(df_long, df_long_Scen3)

df_tot$Metric <- factor(df_tot$Metric, levels = c("MR", "MSEmu1", "MSEmu2", "KLsigma1", "KLsigma2", "MSEmu3", "MSEmu4", "KLsigma3", "KLsigma4"))
metric_labels <- c("MR" = "mMR", 
                   "MSEmu1" = "log[10](MSE*mu[1])",
                   "MSEmu2" = "log[10](MSE*mu[2])",
                   "MSEmu3" = "log[10](MSE*mu[3])",
                   "MSEmu4" = "log[10](MSE*mu[4])",
                   "KLsigma1" = "log[10](KL*Sigma[1])",
                   "KLsigma2" = "log[10](KL*Sigma[2])",
                   "KLsigma3" = "log[10](KL*Sigma[3])",
                   "KLsigma4" = "log[10](KL*Sigma[4])")

# MR
df_tot_MR <- df_tot[which(df_tot$Metric == "MR"), ]
p <- ggplot(df_tot_MR, aes(x = Outlier, y = Value, color = Model, group = Model, shape = Model, linetype = Model)) +
  geom_line() + 
  geom_point() +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 7, 8, 9)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1234", "F2")) +
  facet_wrap(vars(Scenario), scales = "free",
          labeller = labeller(Metric = as_labeller(label_parsed))) +  
  scale_x_discrete(breaks = c(0, 5, 10), labels = c(0, 5, 10)) +  
  xlab("% of outlying values") +
  ylab("mMR") +
  theme(strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5))
p

# PARAMETERS
df_tot_param <- df_tot[which(df_tot$Metric != "MR"), ]
df_tot_param <- df_tot_param %>%
  mutate(ScenarioMetric = paste(Scenario, Metric, sep = "_"))
df_tot_param$ScenarioMetric <- factor(df_tot_param$ScenarioMetric, levels = c(
  "Scenario 1_MSEmu1", "Scenario 1_MSEmu2", "Scenario 2_MSEmu1", "Scenario 2_MSEmu2",
  "Scenario 3_MSEmu1", "Scenario 3_MSEmu2", "Scenario 3_MSEmu3", "Scenario 3_MSEmu4",
  "Scenario 1_KLsigma1", "Scenario 1_KLsigma2", "Scenario 2_KLsigma1", "Scenario 2_KLsigma2",
  "Scenario 3_KLsigma1", "Scenario 3_KLsigma2", "Scenario 3_KLsigma3", "Scenario 3_KLsigma4"))

df_tot_param <- df_tot_param %>%
    mutate(Type = ifelse(grepl("^MSE", Metric), "MSE", "KL"))
df_MSE <- df_tot_param %>% filter(Type == "MSE")
df_KL  <- df_tot_param %>% filter(Type == "KL")

custom_labels <- c(
  "Scenario 1_MSEmu1" = "Scenario 1 \n g = 1",
  "Scenario 1_MSEmu2" = "Scenario 1 \n g = 2",
  "Scenario 2_MSEmu1" = "Scenario 2 \n g = 1",
  "Scenario 2_MSEmu2" = "Scenario 2 \n g = 2",
  "Scenario 3_MSEmu1" = "Scenario 3 \n g = 1",
  "Scenario 3_MSEmu2" = "Scenario 3 \n g = 2",
  "Scenario 3_MSEmu3" = "Scenario 3 \n g = 3",
  "Scenario 3_MSEmu4" = "Scenario 3 \n g = 4",
  "Scenario 1_KLsigma1" = "Scenario 1 \n g = 1",
  "Scenario 1_KLsigma2" = "Scenario 1 \n g = 2",
  "Scenario 2_KLsigma1" = "Scenario 2 \n g = 1",
  "Scenario 2_KLsigma2" = "Scenario 2 \n g = 2",
  "Scenario 3_KLsigma1" = "Scenario 3 \n g = 1",
  "Scenario 3_KLsigma2" = "Scenario 3 \n g = 2",
  "Scenario 3_KLsigma3" = "Scenario 3 \n g = 3",
  "Scenario 3_KLsigma4" = "Scenario 3 \n g = 4"
)

df_MSE_row1 <- df_MSE %>%
  filter(ScenarioMetric %in% c("Scenario 1_MSEmu1", "Scenario 1_MSEmu2", 
                               "Scenario 2_MSEmu1", "Scenario 2_MSEmu2"))
pMSE_row1 <- ggplot(df_MSE_row1, aes(x = Outlier, y = Value, color = Model, group = Model, shape = Model, linetype = Model)) +
  geom_line() + 
  geom_point() +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 7, 8, 9)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1234", "F2")) +
  facet_wrap(vars(ScenarioMetric), scales = "free", nrow = 1, 
             labeller = labeller(ScenarioMetric = as_labeller(custom_labels))) +  
  scale_y_continuous(trans = "log10", labels = scales::label_number()) +
  labs(x = "", y = expression(MSE*mu[g])) +
  theme(strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "none")
df_MSE_row2 <- df_MSE %>%
  filter(!ScenarioMetric %in% c("Scenario 1_MSEmu1", "Scenario 1_MSEmu2", 
                                "Scenario 2_MSEmu1", "Scenario 2_MSEmu2"))
pMSE_row2 <- ggplot(df_MSE_row2, aes(x = Outlier, y = Value, color = Model, group = Model, shape = Model, linetype = Model)) +
  geom_line() + 
  geom_point() +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 7, 8, 9)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1234", "F2")) +
  facet_wrap(vars(ScenarioMetric), scales = "free", nrow = 1, 
             labeller = labeller(ScenarioMetric = as_labeller(custom_labels))) +  
  scale_y_continuous(trans = "log10", labels = scales::label_number()) +
  labs(x = "", y = expression(MSE*mu[g]),) +
  theme(strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "none")

df_KL_row1 <- df_KL %>%
  filter(ScenarioMetric %in% c("Scenario 1_KLsigma1", "Scenario 1_KLsigma2",
                               "Scenario 2_KLsigma1", "Scenario 2_KLsigma2"))
pKL_row1 <- ggplot(df_KL_row1, aes(x = Outlier, y = Value, color = Model, group = Model, shape = Model, linetype = Model)) +
  geom_line() + 
  geom_point() +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 7, 8, 9)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1234", "F2")) +
  facet_wrap(vars(ScenarioMetric), scales = "free", nrow = 1, 
             labeller = labeller(ScenarioMetric = as_labeller(custom_labels))) +
  scale_y_continuous(trans = "log10", labels = scales::label_number()) +
  labs(x = "", y = expression(KL*Sigma[g])) +
  theme(strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "none") 

df_KL_row2 <- df_KL %>%
  filter(!ScenarioMetric %in% c("Scenario 1_KLsigma1", "Scenario 1_KLsigma2",
                                "Scenario 2_KLsigma1", "Scenario 2_KLsigma2"))
pKL_row2 <- ggplot(df_KL_row2, aes(x = Outlier, y = Value, color = Model, group = Model, shape = Model, linetype = Model)) +
  geom_line() + 
  geom_point() +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 7, 8, 9)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1234", "F2")) +
  facet_wrap(vars(ScenarioMetric), scales = "free", nrow = 1, 
             labeller = labeller(ScenarioMetric = as_labeller(custom_labels))) +
  scale_y_continuous(trans = "log10", labels = scales::label_number()) +
  labs(x = "% of outlying values", y = expression(KL*Sigma[g])) +
  theme(strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "bottom")  

p_combined <- pMSE_row1 / pMSE_row2/ pKL_row1 / pKL_row2
p_combined
