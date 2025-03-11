################################################################################
#        Code for replication - Technometrics Review article TCH-24-160        #
#               Simulation Study - Main Article - Figures 1 and 2              #
################################################################################
## To reproduce an entire figure from the simulation study, as reported 
## in the Main Article or in the Supplementary Material, it is necessary to run 
## all scenarios for each level of contamination or at least one scenario of the 
## simulation study for all three level of contamination (0%, 5%, and 10%).
## We have attach an .RData file containing the results of the three scenarios reported 
## in the Main Article, called "Data for Figure Reproducibility.RData", which can 
## be used to fully reproduce Figures 1 and 2.

## ATTENTION: OPEN THIS SCRIPT FROM THE WORKING DIRECTORY OF "Code for Technometrics Reproducibility".
load("Data for Figure Reproducibility.RData")
stat_ggplot
stat_ggplotscen3

## INSTALL PACKAGES (IF NECESSARY)
# install.packages("ggplot2")
# install.packages("tidyr")
# install.packages("dplyr")
library(ggplot2)
library(tidyr)
library(dplyr)
## PREPARE THE OBJECT FOR ggplot
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
p <- ggplot(df_tot_MR, aes(x = Outlier, y = Value, color = Model, group = Model, shape = Model)) +
  geom_line() + 
  geom_point() +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 7, 8, 9)) +
  facet_wrap(vars(Scenario, Metric), scales = "free", 
             labeller = labeller(Metric = as_labeller(metric_labels, label_parsed))) +  
  scale_x_discrete(breaks = c(0, 5, 10), labels = c(0, 5, 10)) +  
  # scale_y_continuous(trans = "log", labels = scales::label_number())) +
  labs(x = "% of outlying values", y = "") +
  theme(strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5))
p

# PARAMETERS
df_tot_param <- df_tot[which(df_tot$Metric != "MR"), ]
p <- ggplot(df_tot_param, aes(x = Outlier, y = Value, color = Model, group = Model, shape = Model)) +
  geom_line() + 
  geom_point() +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 7, 8, 9)) +
  facet_wrap(vars(Scenario, Metric), scales = "free", 
             labeller = labeller(Metric = as_labeller(metric_labels, label_parsed))) +  
  scale_x_discrete(breaks = c(0, 5, 10), labels = c(0, 5, 10)) +  
  scale_y_continuous(trans = "log10", labels = scales::label_number()) +
  labs(x = "% of outlying values", y = "") +
  theme(strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5))
p







