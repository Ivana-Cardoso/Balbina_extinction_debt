#################################
### EXTINCTION DEBT IN BALBINA ##
#################################
# Ivana Cardoso
# ivanawaters@gmail.com
# Created: 19 August 2026
# Last modified: 21 August 2026
#################################

rm(list = ls())
gc()
options(scipen = 999)
set.seed(13)

library(tidyr)
library(dplyr)
library(iNEXT)
library(ggplot2)
library(nlme)
library(tidyverse)

#### STEP ONE ####
#### If the sampling effort were the same, what would be the number of species 
#### in the continuous forest?
#### Rarefaction curves standardized by the number of individuals
#### Continuous forest in PDBFF and Balbina

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Data")
Aurelio_Silva_comm <- read.csv("comm_Marco.csv", header = TRUE)
Bueno_comm <- read.csv("comm_Bueno.csv", header = TRUE)
Amarante_comm <- read.csv("comm_Amarante.csv", header = TRUE)

# Keep only the continuous forest sites
Aurelio_Silva_comm <- Aurelio_Silva_comm[c(11:12),]
Aurelio_Silva_comm <- bind_cols(
  Aurelio_Silva_comm[, 1:6],
  Aurelio_Silva_comm[, 7:103] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)))

Bueno_comm <- Bueno_comm[c(13:22),]
Bueno_comm <- bind_cols(
  Bueno_comm[, 1:6],
  Bueno_comm[, 7:136] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)))

Amarante_comm <- Amarante_comm[c(13:22),]
Amarante_comm <- bind_cols(
  Amarante_comm[, 1:6],
  Amarante_comm[, 7:120] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)))

Balbina_CF <- bind_rows(Aurelio_Silva_comm, Bueno_comm, Amarante_comm)
Balbina_CF[is.na(Balbina_CF)] <- 0
Balbina_CF <- Balbina_CF[,c(2, 7:114)]
Balbina_CF$site <- as.factor(Balbina_CF$site)

CF_Srichness <- vegan::specnumber(Balbina_CF[,2:109])
Sinitial <- max(CF_Srichness)

col_sp <- colnames(Balbina_CF)[2:109]

Balbina_CF <- Balbina_CF %>%
  group_by(site) %>%
  summarise(
    across(all_of(col_sp), sum, na.rm = TRUE),
    .groups = "drop")

Balbina_CF_all <- Balbina_CF %>%
  summarise(
    site = "Balbina_CF",
    across(-site, sum, na.rm = TRUE))
rm(Amarante_comm, Aurelio_Silva_comm, Bueno_comm, col_sp, Balbina_CF)


Lovejoy_comm <- read.csv("lovejoy.csv")

PDBFF_CF <- Lovejoy_comm %>%
  group_by(species_CBRO_2021) %>%
  summarise(n = sum(band), .groups = "drop") %>%
  pivot_wider(
    names_from = species_CBRO_2021,
    values_from = n,
    values_fill = 0)

PDBFF_CF$site <- "PDBFF_CF"
PDBFF_CF <- PDBFF_CF[,c(144, 1:143)]
rm(Lovejoy_comm)

CFs <- bind_rows(PDBFF_CF, Balbina_CF_all)
CFs[is.na(CFs)] <- 0
rownames(CFs) <- CFs$site

CFs_iNEXT <- list(
  PDBFF = as.numeric(PDBFF_CF[1, -1]),
  Balbina = as.numeric(Balbina_CF_all[1, -1]))

n <- sapply(CFs_iNEXT, sum) # Number of individuals in PDBFF and Balbina
endpoint_min <- 2 * min(n) # Define the endpoint as the double of Balbina

# PDBFF rarefied to Balbina's double of number of individuals
out <- iNEXT(CFs_iNEXT, q = 0, datatype = "abundance", 
             endpoint = endpoint_min,
             nboot = 1000)
table <- out$iNextEst$size_based

# Balbina extrapolated to observed richness in PDBFF
# out2 <- iNEXT(CFs_iNEXT, q = 0, datatype = "abundance", 
#             endpoint = 10916,
#            nboot = 100)
# table2 <- out2$iNextEst$size_based

y <- out$iNextEst$size_based
z <- y[which(y$Method != "Observed"),]

z$Assemblage = factor(z$Assemblage,
                      levels = unique(z$Assemblage))

iNEXT_plot_2xBalbina <-
  ggplot(data = z,
         aes(x = m, y = qD, 
             colour = Assemblage, 
             fill = Assemblage,
             linetype = Method)) +
  labs(y = "Number of species",
       x = "Number of individuals") +
  scale_colour_manual(values = c("PDBFF" = "#D4A72C","Balbina" = "#2C7FB8")) +
  scale_fill_manual(values = c("PDBFF" = "#D4A72C", "Balbina" = "#2C7FB8")) +
  scale_linetype_manual(values = c("Rarefaction" = "solid",
                                   "Extrapolation" = "dashed"))+
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL),
              colour = NA, alpha = 0.30, 
              show.legend = FALSE) +
  geom_line(size = 1.3) +
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +
  scale_x_continuous(breaks = c(seq(0, 2500, by = 500), 2850))+
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = c(0.80, 0.30))
iNEXT_plot_2xBalbina

rm(Balbina_CF_all, CFs, PDBFF_CF, endpoint_min, n, CF_Srichness)
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("iNEXT_Balbina_PDBFF.png",
       plot = iNEXT_plot_2xBalbina,
       width = 18, height = 12,
       units = "cm", dpi = 900)



#### STEP TWO ####
#### Estimating initial species richness (S0) based on PDBFF data: 
#### (Ferraz et al., 2003) Rates of species loss from Amazonian forest
#### fragments. Table 1:
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt")
Balbina <- read.csv("dados_completos_5periodos.csv", row.names = 1)
# Some islands sampled by Marco Aurelio Silva were sampled in 2010, 2011, or both years. This info is in the name of the site. BUT I treated all of them as having the same period (column 3, “periodo_5”), so they had the same isolation time (column 4, “t5”, 23 years) and could be grouped together.

a = c(2.8, 1.6, 1.8, 1.1, 1.6, 13, 11.2, 10.7, 11, 98.1, 101.2)
s = c(90, 73, 86, 84, 82, 85, 101, 92, 89, 115, 111)

m0 = lm(log10(s) ~ log10(a))
summary(m0)

data_debt <- data.frame(a = round(Balbina$area, 2))
colnames(data_debt)[1] <- "a"
rownames(data_debt) <- rownames(Balbina)

data_debt$s_prevista <- 10^(predict(m0, newdata = data.frame(a = data_debt$a)))
data_debt$Sinitial_max <- Sinitial

plot(log10(a), log10(s),
     xlab = "log10(Area)",
     ylab = "log10(Species richness)",
     main = "S0 for Balbina islands")

abline(m0)

points(log10(data_debt$a[seq(1, 82, by = 6)]),
       log10(data_debt$s_prevista[seq(1, 82, by = 6)]),
       pch = 19, cex = 1.5, col = "blue")

text(log10(data_debt$a[seq(1, 82, by = 6)]),
     log10(data_debt$s_prevista[seq(1, 82, by = 6)]),
     labels = paste(round(data_debt$s[seq(1, 82, by = 6)]), "species"),
     pos = 3, col = "blue")

text(log10(data_debt$a[seq(1, 82, by = 6)]),
     log10(data_debt$s_prevista[seq(1, 82, by = 6)]),
     labels = paste(round(data_debt$a[seq(1, 82, by = 6)]), "ha"),
     pos = 1, col = "red")



#### STEP THREE ####
#### Calculating the number of species remaining on each period
#### number of remaining species = Sobserved / Sinitial 

data_debt$site <- rownames(data_debt)
Balbina$site <- rownames(Balbina)

data_debt <- data_debt[,c(4,1:3)]
data_debt <- data_debt %>%
  left_join(Balbina %>% select(site, obs.richness),
            by = "site")

data_debt$"remaining_sp%" <- data_debt$obs.richness/data_debt$Sinitial_max * 100

data_debt <- data_debt %>%
  left_join(Balbina %>% select(site, periodo_5),
            by = "site")
data_debt$periodo_5 <- as.factor(data_debt$periodo_5)

means <- tapply(data_debt$`remaining_sp%`, data_debt$periodo_5, 
                mean, na.rm = TRUE)

remaining_sp_plot <-
  ggplot(mapping = aes(x = periodo_5, y = `remaining_sp%`),
         data = data_debt) +
  labs(x = "Ano de amostragem", y = "Espécies restantes (%)") +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(
    labels = c("2010" = "2010-2011", "2015" = "2015",
               "2016" = "2016", "2023" = "2023",
               "2024" = "2024")) +
  scale_y_continuous(breaks = seq(0, max(data_debt$`remaining_sp%`), by = 10)) +
  geom_point(position = position_jitter(), size = 3, alpha = 0.3) +
  stat_summary(fun = mean, geom = "point",
               size = 4, color = "red") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "bottom")
remaining_sp_plot

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("remaining_sp.png",
       plot = remaining_sp_plot,
       width = 14, height = 12,
       units = "cm", dpi = 900)



#### STEP FOUR ####
#### Fit extinction debt model 

Balbina_m1 <- subset(Balbina, obs.richness > 0) #two islands has 0 species
m1 = lm(log10(obs.richness) ~ log10(area), data = Balbina_m1)
summary(m1)

c0 <- 10^(coef(m0)[1]) #check the c value from before isolation
z0 <- coef(m0)[2] #check the z value from before isolation

cc <- 10^(coef(m1)[1]) #check the c value of all 5 periods sampled
zz <- coef(m1)[2] #check the z value of all 5 periods samples

colnames(data_debt)[2] <- "A"
colnames(data_debt)[3] <- "S0"
data_debt$periodo_5 <- as.numeric(as.character(data_debt$periodo_5))
data_debt$t <- data_debt$periodo_5 - 1987
#data_debt <- subset(data_debt, obs.richness > 0)
data_debt$ilhas <- sub("_.*", "", data_debt$site)


###############################
#### NEGATIVE EXPONENTIAL MODEL
#### S(t,A) = cAz + (S0 - cAz)exp(-t / (k0 * A^alpha))
####
#### S(t,A) is observed richness on time t
#### cAz is the number of species in equilibrium
#### S0 is the number of species before isolation - max of Balbina's CFs
#### t is the time since isolation
#### k0 is the constant of relaxation time
#### alpha is the scaling exponent for the relaxation time with respect to area; it indicates how the speed at which the extinction debt is paid off changes as the island varies in size. Halley uses alpha ~ 0.5

debt_ne <- nlme(
  obs.richness ~ c * A^z + (Sinitial_max - c * A^z) * exp(-t / (k0 * A^alpha)),
  data = data_debt,
  fixed = c + z + k0 + alpha ~ 1,
  random = k0 ~ 1 | ilhas,
  correlation = corCAR1(form = ~ t | ilhas),
  start = c(cc, zz, 10, 0.5),
  control = nlmeControl(maxIter = 500, msMaxIter = 500))

summary(debt_ne)
intervals(debt_ne, which = "fixed")



###############################
#### POWER LAW MODEL
#### S(t,A) = cAz + (S0 - cAz)*(1 + t / (k0 * A^alpha))^(-beta)
####
#### S(t,A) is observed richness on time t
#### cAz is the number of species in equilibrium
#### S0 is the number of species before isolation - max of Balbina's CFs
#### t is the time since isolation
#### k0 is the constant of relaxation time
#### alpha is the scaling exponent for the relaxation time with respect to area; it indicates how the speed at which the extinction debt is paid off changes as the island varies in size. Halley uses alpha ~ 0.5
#### beta controls the shape of the decay tail; higher values indicate faster relaxation and shorter tails, whereas lower values indicate slower, more prolonged relaxation.

debt_pl <- nlme(
  obs.richness ~ c * A^z + (Sinitial_max - c * A^z) *
    (1 + t / (k0 * A^alpha))^(-beta),
  data = data_debt,
  fixed = c + z + k0 + alpha + beta ~ 1,
  random = k0 ~ 1 | ilhas,
  correlation = corCAR1(form = ~ t | ilhas),
  start = c(cc, zz, 10, 0.5, 0.8),
  control = nlmeControl(maxIter = 500,msMaxIter = 500))

summary(debt_pl)
intervals(debt_pl)
AIC(debt_ne, debt_pl)  #negative exp is better for the data (delta = 2.79)


#### I will keep using the negative exponential model below

# Calculating pseudo-R2
y_obs <- fitted(debt_ne) + residuals(debt_ne)
R2 <- 1 - sum(residuals(debt_ne)^2) /
  sum((y_obs - mean(y_obs))^2)
R2 # 0.67

# Extracting the coefficients from the model - fixed effects
fixed_coefs <- fixef(debt_ne)
c_coef     <- fixed_coefs["c"]
z_coef     <- fixed_coefs["z"]
alpha_coef <- fixed_coefs["alpha"]

# Island-specific k0 (fixed effect + island-level random effect)
coef_por_ilha <- coef(debt_ne) %>%
  rownames_to_column("ilhas") %>%
  select(ilhas, k0)   

# Extracting island-specific parameters
ilha_params <- data_debt %>%
  group_by(ilhas) %>%
  summarise(
    A = first(A),
    Sinitial_mean = first(Sinitial_max),
    .groups = "drop"
  ) %>%
  left_join(coef_por_ilha, by = "ilhas") %>%
  rename(k0_ilha = k0) %>%
  mutate(
    k_ilha = k0_ilha * A^alpha_coef,
    t_half_life = k_ilha * log(2))

# Mean half-life between islands
cat("Mean half-life between islands:",
    mean(ilha_params$t_half_life), "\n")

t_half_ilhas <- ilha_params %>%
  select(ilhas, A, k0_ilha, k_ilha, t_half_life)
print(t_half_ilhas)

# Calculating the extinction debt
data_debt$S_eq <- c_coef * data_debt$A^z_coef
data_debt$divida_inicial <-
  data_debt$Sinitial_max - data_debt$S_eq

data_debt <- data_debt %>%
  left_join(ilha_params %>% select(ilhas, k_ilha),
            by = "ilhas")

data_debt$divida_restante <-
  data_debt$divida_inicial *
  exp(-data_debt$t / data_debt$k_ilha)

# Relaxation curves by island
t_seq <- seq(0, max(data_debt$t) * 2, length.out = 200)

curvas <- ilha_params %>%
  crossing(t = t_seq) %>%
  mutate(S_eq = c_coef * A^z_coef,
         S_pred = S_eq + (Sinitial_mean - S_eq) * exp(-t / k_ilha))

biogeograp_model_plot <-
  ggplot() +
  geom_line(data = curvas,
            aes(x = t, y = S_pred, group = ilhas, colour = A),
            linewidth = 1) +
  
  geom_point(data = data_debt, aes(x = t, y = obs.richness),
             colour = "black", size = 2) +
  
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  
  scale_colour_viridis_c(name = "Island area (ha)") +
  
  labs(x = "Time since isolation (years)",
       y = "Number of species remaining") +
  
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "right")

biogeograp_model_plot

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("biogeograp_model_plot.png",
       plot = biogeograp_model_plot,
       width = 24,
       height = 12,
       units = "cm",
       dpi = 900)




#### VISUALIZING DATA ####
# SARs for S0 derived from PDBFF data, and for the five sampled periods in Balbina
models_plot <- 
  ggplot() +
  geom_smooth(aes(x = log10(a), y = log10(s), color = "SAR para S0"),
              method = "lm", se = TRUE, fill = "#D4A72C") +
  geom_point(aes(x = log10(a), y = log10(s), color = "SAR para S0"),
             size = 3, alpha = 0.4) +
  geom_smooth(data = Balbina_m1, 
              aes(x = log10(area), y = log10(obs.richness), 
                  color = "SAR para St"),
              method = "lm", se = TRUE, fill = "#2C7FB8") +
  geom_point(data = Balbina_m1, aes(x = log10(area), y = log10(obs.richness), 
                                    color = "SAR para St"), 
             size = 3, alpha = 0.4) +
  scale_color_manual(name = NULL, values = c("SAR para S0" = "#D4A72C",
                                             "SAR para St" = "#2C7FB8")) +
  xlim(0, 3) +
  ylim(0, 3) +
  labs(x = "Area (log10)", y = "Number of species (log10)") +
  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "bottom")
models_plot

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("models_S0_st.png",
       plot = models_plot,
       width = 12,
       height = 12,
       units = "cm",
       dpi = 900)

#Initial, Observed and in Equilibrium Richness
plot_data <- data_debt %>%
  select(ilhas, Sinitial_max, obs.richness, S_eq) %>%
  pivot_longer(
    cols = c(Sinitial_max, obs.richness, S_eq),
    names_to = "stage",
    values_to = "richness"
  ) %>%
  mutate(
    stage = factor(
      stage,
      levels = c("Sinitial_max", "obs.richness", "S_eq"),
      labels = c("Initial", "Observed", "Equilibrium")
    )
  )

Richness_plot <-
ggplot(plot_data, aes(x = stage, y = richness, group = ilhas)) +
  geom_line(linewidth = 0.2) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = NULL, y = "Species richness") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "none")
Richness_plot

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("initial_obs_equilibrium_richness.png",
       plot = Richness_plot,
       width = 12,
       height = 8,
       units = "cm",
       dpi = 900)