#################################
### EXTINCTION DEBT IN BALBINA ##
#################################
# Ivana Cardoso
# ivanawaters@gmail.com
# Created: 19 August 2026
# Last modified: 25 August 2026
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

#########################
#### ORGANIZING DATA ####
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Data")
Aurelio_Silva_comm <- read.csv("comm_Marco.csv", header = TRUE)
Bueno_comm <- read.csv("comm_Bueno.csv", header = TRUE)
Amarante_comm <- read.csv("comm_Amarante.csv", header = TRUE)

Aurelio_Siva_env <- read.csv("environment_data_Marco2010.csv", header = TRUE)
Bueno_env <- read.csv("environment_data_Anderson2015.csv", header = TRUE)
Amarante_env <- read.csv("environment_data_Amarante2023.csv", header = TRUE)

# Keep only the continuous forest sites
Aurelio_Silva_comm_CF <- Aurelio_Silva_comm[c(12:14),]
Aurelio_Silva_comm_CF <- bind_cols(
  Aurelio_Silva_comm_CF[, 1:7],
  Aurelio_Silva_comm_CF[, 8:104] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)))

Bueno_comm_CF <- Bueno_comm[c(13:22),]
Bueno_comm_CF <- bind_cols(
  Bueno_comm_CF[, 1:7],
  Bueno_comm_CF[, 8:136] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)))

Amarante_comm_CF <- Amarante_comm[c(13:22),]
Amarante_comm_CF <- bind_cols(
  Amarante_comm_CF[, 1:7],
  Amarante_comm_CF[, 8:120] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)))

mean(c(
  vegan::specnumber(Aurelio_Silva_comm_CF[,-c(1:7)]),
  vegan::specnumber(Bueno_comm_CF[,-c(1:6)]),
  vegan::specnumber(Amarante_comm_CF[,-c(1:6)])))

Balbina_CF <- bind_rows(Aurelio_Silva_comm_CF, Bueno_comm_CF, Amarante_comm_CF)
Balbina_CF[is.na(Balbina_CF)] <- 0
Balbina_CF <- Balbina_CF[,c(2, 8:115)]
Balbina_CF$Local <- as.factor(Balbina_CF$Local)

colnames(Balbina_CF)[1] <- "site"
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

# Keep only the islands
Aurelio_Silva_comm_islands <- Aurelio_Silva_comm[-c(12:14),]
Aurelio_Silva_comm_islands <- bind_cols(
  Aurelio_Silva_comm_islands[, 1:7],
  Aurelio_Silva_comm_islands[, 8:104] %>% 
    select(where(~ sum(.x, na.rm = TRUE) > 0)))
Aurelio_Silva_comm_islands <- Aurelio_Silva_comm_islands %>%
  left_join(Aurelio_Siva_env %>% select(Local, area),
            by = "Local")
Aurelio_Silva_comm_islands <- Aurelio_Silva_comm_islands[,c(1:3,85,4:84)]
Aurelio_Silva_comm_islands$area[c(10,11)] <- Aurelio_Siva_env$area[25]
Aurelio_Silva_comm_islands$area[c(14,15)] <- Aurelio_Siva_env$area[10]
Aurelio_Silva_comm_islands$area[c(16,17)] <- Aurelio_Siva_env$area[22]

# mean transects in the same island
plots_mean <- c("Beco_do_catitu", "Gaviao_real", "Martelo")

island_mean <- Aurelio_Silva_comm_islands %>%
  filter(ilha %in% plots_mean) %>%
  group_by(ilha) %>%
  summarise(across(1:8, first),
            across(9:last_col(), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")
island_mean[,c(9:85)] <- round(island_mean[,c(9:85)],0)

# summing islands that were sampled two years
parte_resto <- Aurelio_Silva_comm_islands %>%
  filter(!ilha %in% plots_mean) %>%
  group_by(ilha) %>%
  summarise(across(1:8, first),
            across(9:last_col(), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

Aurelio_Silva_comm_islands <- bind_rows(island_mean, parte_resto)
Aurelio_Silva_comm_islands$ilha[1] <- "Beco_do_Catitu"
Aurelio_Silva_comm_islands$Ilha_Ano[1] <- "Beco_do_Catitu_B_2010"
Aurelio_Silva_comm_islands$Local[1] <- "Beco_do_Catitu_B"

Bueno_comm_islands <- Bueno_comm[-c(13:22),]
Bueno_comm_islands <- bind_cols(
  Bueno_comm_islands[, 1:6],
  Bueno_comm_islands[, 7:136] %>% 
    select(where(~ sum(.x, na.rm = TRUE) > 0)))
colnames(Bueno_env)[1] <- "Local"
Bueno_comm_islands <- Bueno_comm_islands %>%
  left_join(Bueno_env %>% select(Local, area),
            by = "Local")
Bueno_comm_islands <- Bueno_comm_islands[,c(1:2,116,3:115)]

Amarante_comm_islands <- Amarante_comm[-c(13:22),]
Amarante_comm_islands <- bind_cols(
  Amarante_comm_islands[, 1:6],
  Amarante_comm_islands[, 7:120] %>% 
    select(where(~ sum(.x, na.rm = TRUE) > 0)))
colnames(Amarante_env)[1] <- "Local"
Amarante_comm_islands <- Amarante_comm_islands %>%
  left_join(Amarante_env %>% select(Local, area),
            by = "Local")
Amarante_comm_islands <- Amarante_comm_islands[,c(1:2,96,3:95)]

Balbina_islands <- bind_rows(Aurelio_Silva_comm_islands, Bueno_comm_islands, Amarante_comm_islands)
Balbina_islands$ilha[20:149] <- Balbina_islands$Local[20:149]
Balbina_islands <- as.data.frame(Balbina_islands)
Balbina_islands[is.na(Balbina_islands)] <- 0

Balbina_islands$obs.richness <- vegan::specnumber(Balbina_islands[,9:151])

##### ONLY THE 11 ISLANDS SAMPLED DURING THE 5 PERIODS
ilhas_Marco <- Balbina_islands %>%
  filter(Ano %in% c(2010, 2011)) %>%
  distinct(ilha) %>%
  pull(ilha)

Balbina_11_islands <- Balbina_islands %>%
  filter(ilha %in% ilhas_Marco)

Balbina_11_islands <- Balbina_11_islands[-c(1:8),]

rm(Aurelio_Silva_comm_islands, Aurelio_Silva_comm, Bueno_comm, 
   Bueno_comm_islands, Amarante_comm_islands, Amarante_comm,
   Aurelio_Siva_env, Bueno_env, Amarante_env, Aurelio_Silva_comm_CF,
   Bueno_comm_CF, Amarante_comm_CF)



#### STEP ONE ####
#### If the sampling effort were the same, what would be the number of species 
#### in the continuous forest?
#### Rarefaction curves standardized by the number of individuals
#### Continuous forest in PDBFF and Balbina
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

rm(Balbina_CF_all, CFs, PDBFF_CF, endpoint_min, n, CF_Srichness, Lovejoy_comm)
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("iNEXT_Balbina_PDBFF.png",
       plot = iNEXT_plot_2xBalbina,
       width = 18, height = 12,
       units = "cm", dpi = 900)


#### STEP TWO ####
#### Estimating initial species richness (S0) based on PDBFF data: 
#### (Ferraz et al., 2003) Rates of species loss from Amazonian forest
#### fragments. Table 1:

a = c(2.8, 1.6, 1.8, 1.1, 1.6, 13, 11.2, 10.7, 11, 98.1, 101.2)
s = c(90, 73, 86, 84, 82, 85, 101, 92, 89, 115, 111)

m0 = lm(log10(s) ~ log10(a))
summary(m0)

data_debt <- data.frame(a = round(Balbina_islands$area, 2))
colnames(data_debt)[1] <- "a"
rownames(data_debt) <- Balbina_islands$Ilha_Ano

# using the SAR from PDBFF sampled before isolation to predict what would be the species richnes in Balbina, in areas with the same size of the islands, when it was not isolated yet
data_debt$s_prevista <- 10^(predict(m0, newdata = data.frame(a = data_debt$a))) 

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
#### Calculating the number of species on each period

data_debt$site <- rownames(data_debt)
Balbina_islands$site <- Balbina_islands$Ilha_Ano

data_debt <- data_debt[,c(3, 1:2)]
data_debt <- data_debt %>%
  left_join(Balbina_islands %>% select(site, obs.richness),
            by = "site")
data_debt <- data_debt %>%
  left_join(Balbina_islands %>% select(site, t),
            by = "site")
data_debt <- data_debt %>%
  left_join(Balbina_islands %>% select(site, t_period),
            by = "site")
data_debt <- data_debt %>%
  left_join(Balbina_islands %>% select(site, Ano),
            by = "site")
data_debt$Ano[1:19] <- 2010 #Im considering Marco's as one period
data_debt$Ano <- as.factor(data_debt$Ano)

means <- tapply(data_debt$obs.richness, data_debt$Ano, 
                mean, na.rm = TRUE)
means

years_sp_plot <-
  ggplot(mapping = aes(x = Ano, y = obs.richness),
         data = data_debt) +
  labs(x = "Sampling year", y = "Number of species") +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(
    labels = c("2010" = "2010-2011", "2015" = "2015",
               "2016" = "2016", "2023" = "2023",
               "2024" = "2024")) +
  scale_y_continuous(breaks = seq(0, max(data_debt$obs.richness), by = 10)) +
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
years_sp_plot

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("years_sp.png",
       plot = years_sp_plot,
       width = 14, height = 12,
       units = "cm", dpi = 900)


#### STEP FOUR ####
#### Fit extinction debt model 

Balbina_m1 <- subset(data_debt, data_debt$obs.richness > 0) #two islands has 0 species
m1 = lm(log10(obs.richness) ~ log10(a), data = Balbina_m1)
summary(m1)

cc <- 10^(coef(m1)[1]) #check the c value of all 5 periods sampled
zz <- coef(m1)[2] #check the z value of all 5 periods samples

colnames(data_debt)[2] <- "A"
data_debt$ilhas <- sub("(_[A-Za-z])?_[0-9]{4}$", "", data_debt$site)

###############################
#### NEGATIVE EXPONENTIAL MODEL
#### S(t,A) = cAz + (S0 - cAz)exp(-t / (k0 * A^alpha))
####
#### S(t,A) is observed richness on time t
#### cAz is the number of species in equilibrium
#### S0 is the number of species before isolation
#### t is the time since isolation
#### k0 is the constant of relaxation time
#### alpha is the scaling exponent for the relaxation time with respect to area; it indicates how the speed at which the extinction debt is paid off changes as the island varies in size. Halley et al., 2016 uses alpha ~ 0.5
debt_ne <- nlme(
  obs.richness ~ c * A^z + (s_prevista - c * A^z) * exp(-t / (k0 * A^alpha)),
  data = data_debt,
  fixed = c + z + k0 + alpha ~ 1,
  random = k0 ~ 1 | ilhas,
  correlation = corCAR1(form = ~ t | ilhas),
  start = c(cc, zz, 10, 0.52),
  control = nlmeControl(maxIter = 500, msMaxIter = 500))

summary(debt_ne)
intervals(debt_ne, which = "fixed")


#### Creating the functions based on the model ####
#### Below is a function, following the formula, to calculate de number of species predicted to occur on time t and area A (I need to specify). I also need to specify the S0, but I will derivate from PDBFF
calc_S0 <- function(A) {
  10^predict(m0, newdata = data.frame(a = A))}

coefs <- fixef(debt_ne)
S_pred <- function(A, t, S0,
                   c = coefs["c"],
                   z = coefs["z"],
                   k0 = coefs["k0"],
                   alpha = coefs["alpha"]) {
  c * A^z + (S0 - c * A^z) * exp(-t / (k0 * A^alpha))
}

#### lets test it with 12.64 ha
S_pred(A = 12.64, t = 39, S0 = calc_S0(A = 12.64))


#### Visualizing the model for the 41 islands ####
# Calculating pseudo-R2
y_obs <- fitted(debt_ne) + residuals(debt_ne)
R2 <- 1 - sum(residuals(debt_ne)^2) /
  sum((y_obs - mean(y_obs))^2)
R2 # 0.60

# Extracting the coefficients from the model - fixed effects
c_coef <- coefs["c"]
z_coef <- coefs["z"]
k0_coef <- coefs["k0"]
alpha_coef <- coefs["alpha"]

# Island-specific k0 (fixed effect + island-level random effect)
coef_por_ilha <- coef(debt_ne) %>%
  rownames_to_column("ilhas") %>%
  select(ilhas, k0)   

# Extracting island-specific parameters
ilha_params <- data_debt %>%
  group_by(ilhas) %>%
  summarise(
    A = last(A),
    s_prevista = mean(s_prevista),
    .groups = "drop"
  ) %>%
  left_join(coef_por_ilha, by = "ilhas") %>%
  rename(k0_ilha = k0) %>%
  mutate(
    k_ilha = k0_ilha * A^alpha_coef,
    t_half_life = k_ilha * log(2))


##### Métrica = Extinction debt paid (%)
##### debt_paid = (S_0 - S_2026) / (S_0 - S_eq) * 100
richness_2026 <- ilha_params %>%
  mutate(S_eq = c_coef * A^z_coef, 
         S_pred_2026 = S_eq + 
           (s_prevista - S_eq) * exp(-39 / k_ilha)) %>%
  select(ilhas, A, s_prevista, k_ilha, S_eq, S_pred_2026)
print(richness_2026)

debt_paid <- (richness_2026$s_prevista - richness_2026$S_pred_2026) / 
  (richness_2026$s_prevista - richness_2026$S_eq) * 100

paste0("Estimamos que uma ilha de ", richness_2026$A, 
       " ha já tenha pago ", round(debt_paid), 
       "% de seu débito de extinção após ", 2026 - 1987, 
       " anos de isolamento.")

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
png("debt_paid_area.png", width = 14, height = 10,units = "cm", res = 900)
plot(debt_paid ~ log10(richness_2026$A),
     main = "Paid extinction debt over 39 years of island creation",
     xlab = "Island area (ha)",
     ylab = "Extinction debt paid (%)",
     xaxt = "n",
     las = 1,
     pch = 21, col = "black", bg = "lightgrey", cex = 1.5)+
  axis(1,
       at = log10(c(1, 10, 100, 1000)),
       labels = c(1, 10, 100, 1000))+
  abline(v = log10(c(1, 10, 100, 1000)),
         lty = "dashed", col = "lightgrey")
dev.off()


# Mean half-life between islands
cat("Mean half-life between islands:",
    mean(ilha_params$t_half_life), "\n")
max(ilha_params$t_half_life)
min(ilha_params$t_half_life)

t_half_ilhas <- ilha_params %>%
  select(ilhas, A, k0_ilha, k_ilha, t_half_life)
print(t_half_ilhas)

estimar_meia_vida <- function(A) {
  k <- k0_coef * A^alpha_coef
  t_half <- k * log(2)
  data.frame(A = A, k = k, t_half_life = t_half)}

estimar_meia_vida(100)


# Relaxation curves by island
t_seq <- seq(0, 100, length.out = 200)

curvas <- ilha_params %>%
  crossing(t = t_seq) %>%
  mutate(S_eq = c_coef * A^z_coef,
         S_pred = S_eq + (s_prevista - S_eq) * exp(-t / k_ilha))

biogeograp_model_plot <-
  ggplot() +
  geom_line(data = curvas,
            aes(x = t, y = S_pred, group = ilhas, colour = A),
            linewidth = 1) +
   geom_point(
    data = data_debt,
    aes(x = t, y = obs.richness, colour = A),
    size = 2, alpha = 0.6,
    position = position_jitter(width = 2, height = 0)) +
  
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  
  scale_colour_viridis_c(name = "Island area (ha)") +
  
  labs(x = "Time since isolation (years)",
       y = "Number of species remaining") +
  
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size = 10))
biogeograp_model_plot

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("biogeograp_model.png",
       plot = biogeograp_model_plot,
       width = 24,
       height = 12,
       units = "cm",
       dpi = 900)


####################################
##### Plot with reference areas ####

# Indicating the reference values for area
areas_ref <- c(1, 10, 100, 1000)
s_prevista_ref <- 10^(predict(m0, newdata = data.frame(a = areas_ref)))

# Predicted curves
curvas_ref <- expand.grid(A = areas_ref, t = t_seq) %>%
  mutate(
    S_eq   = c_coef * A^z_coef,
    k      = k0_coef * A^alpha_coef,
    S_pred = S_eq + (s_prevista_ref - S_eq) * exp(-t / k),
    A = factor(A, levels = areas_ref))

biogeograp_model_plot_ref <-
  ggplot() +
  geom_line(data = curvas,
            aes(x = t, y = S_pred, group = ilhas),
            colour = "grey85", linewidth = 0.4) +
  # geom_point(data = data_debt, aes(x = t, y = obs.richness),
  #           colour = "black", size = 2, alpha = 0.4) +
  geom_line(data = curvas_ref,
            aes(x = t, y = S_pred, group = A, colour = A),
            linewidth = 1) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  scale_colour_viridis_d(name = "Island area (ha)") +
  labs(x = "Time since isolation (years)",
       y = "Number of species") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size = 10))

biogeograp_model_plot_ref
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("biogeograp_model_reference.png",
       plot = biogeograp_model_plot_ref,
       width = 24,
       height = 12,
       units = "cm",
       dpi = 900)


#### FOR THE 11 ISLANDS SAMPLED OVER 5 PERIODS ####
#### Estimating initial species richness (S0) based on PDBFF data: 
a = c(2.8, 1.6, 1.8, 1.1, 1.6, 13, 11.2, 10.7, 11, 98.1, 101.2)
s = c(90, 73, 86, 84, 82, 85, 101, 92, 89, 115, 111)

m0 = lm(log10(s) ~ log10(a))
summary(m0)

data_debt_11 <- data.frame(a = round(Balbina_11_islands$area, 2))
colnames(data_debt_11)[1] <- "a"
rownames(data_debt_11) <- Balbina_11_islands$Ilha_Ano

data_debt_11$s_prevista <- 10^(predict(m0, newdata = data.frame(a = data_debt_11$a))) 

#### Calculating the number of species on each period
data_debt_11$site <- rownames(data_debt_11)
Balbina_11_islands$site <- Balbina_11_islands$Ilha_Ano

data_debt_11 <- data_debt_11[,c(3, 1:2)]
data_debt_11 <- data_debt_11 %>%
  left_join(Balbina_11_islands %>% select(site, obs.richness),
            by = "site")
data_debt_11 <- data_debt_11 %>%
  left_join(Balbina_11_islands %>% select(site, t),
            by = "site")
data_debt_11 <- data_debt_11 %>%
  left_join(Balbina_11_islands %>% select(site, t_period),
            by = "site")
data_debt_11 <- data_debt_11 %>%
  left_join(Balbina_11_islands %>% select(site, Ano),
            by = "site")
data_debt_11$Ano[1:11] <- 2010 #Im considering Marco's as one period
data_debt_11$Ano <- as.factor(data_debt_11$Ano)

means <- tapply(data_debt_11$obs.richness, data_debt_11$Ano, 
                mean, na.rm = TRUE)
means

years_sp_plot <-
  ggplot(mapping = aes(x = Ano, y = obs.richness),
         data = data_debt_11) +
  labs(x = "Sampling year", y = "Number of species") +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(
    labels = c("2010" = "2010-2011", "2015" = "2015",
               "2016" = "2016", "2023" = "2023",
               "2024" = "2024")) +
  scale_y_continuous(breaks = seq(0, max(data_debt_11$obs.richness), by = 10)) +
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
years_sp_plot


#### Fit extinction debt model 
rownames(Balbina_11_islands) <- Balbina_11_islands$Ilha_Ano
data_11 <- as.data.frame(t(Balbina_11_islands[,c(9:151)]))
data_11 <- Filter(function(x) sum(x) != 0, data_11)

Balbina_m1 <- subset(data_debt_11, data_debt_11$obs.richness > 0) #two islands has 0 species
m1 = lm(log10(obs.richness) ~ log10(a), data = Balbina_m1)
summary(m1)

cc <- 10^(coef(m1)[1]) #check the c value of all 5 periods sampled
zz <- coef(m1)[2] #check the z value of all 5 periods samples

colnames(data_debt_11)[2] <- "A"
data_debt_11$ilhas <- sub("(_[A-Za-z])?_[0-9]{4}$", "", data_debt_11$site)

#### NEGATIVE EXPONENTIAL MODEL

debt_ne <- nlme(
  obs.richness ~ c * A^z + (s_prevista - c * A^z) * exp(-t / (k0 * A^alpha)),
  data = data_debt_11,
  fixed = c + z + k0 + alpha ~ 1,
  random = k0 ~ 1 | ilhas,
  correlation = corCAR1(form = ~ t | ilhas),
  start = c(cc, zz, 5, 0.52),
  control = nlmeControl(maxIter = 500, msMaxIter = 500))

summary(debt_ne)
intervals(debt_ne, which = "fixed")

#### Creating the functions based on the model 
calc_S0 <- function(A) {
  10^predict(m0, newdata = data.frame(a = A))}

coefs <- fixef(debt_ne)
S_pred <- function(A, t, S0,
                   c = coefs["c"],
                   z = coefs["z"],
                   k0 = coefs["k0"],
                   alpha = coefs["alpha"]) {
  c * A^z + (S0 - c * A^z) * exp(-t / (k0 * A^alpha))
}

#### lets test it with 100 ha
S_pred(A = 100, t = 40, S0 = calc_S0(A = 100))


#### Visualizing the model for the 11 islands ####
# Calculating pseudo-R2
y_obs <- fitted(debt_ne) + residuals(debt_ne)
R2 <- 1 - sum(residuals(debt_ne)^2) /
  sum((y_obs - mean(y_obs))^2)
R2 # 0.45

# Extracting the coefficients from the model - fixed effects
c_coef <- coefs["c"]
z_coef <- coefs["z"]
k0_coef <- coefs["k0"]
alpha_coef <- coefs["alpha"]

# Island-specific k0 (fixed effect + island-level random effect)
coef_por_ilha <- coef(debt_ne) %>%
  rownames_to_column("ilhas") %>%
  select(ilhas, k0)   

# Extracting island-specific parameters
ilha_params <- data_debt_11 %>%
  group_by(ilhas) %>%
  summarise(
    A = last(A),
    s_prevista = mean(s_prevista),
    .groups = "drop"
  ) %>%
  left_join(coef_por_ilha, by = "ilhas") %>%
  rename(k0_ilha = k0) %>%
  mutate(
    k_ilha = k0_ilha * A^alpha_coef,
    t_half_life = k_ilha * log(2))


##### Métrica = Extinction debt paid (%)
##### debt_paid = (S_0 - S_2026) / (S_0 - S_eq) * 100
richness_2026 <- ilha_params %>%
  mutate(S_eq = c_coef * A^z_coef, 
         S_pred_2026 = S_eq + 
           (s_prevista - S_eq) * exp(-39 / k_ilha)) %>%
  select(ilhas, A, s_prevista, k_ilha, S_eq, S_pred_2026)
print(richness_2026)

debt_paid <- (richness_2026$s_prevista - richness_2026$S_pred_2026) / 
  (richness_2026$s_prevista - richness_2026$S_eq) * 100

paste0("Estimamos que uma ilha de ", richness_2026$A, 
       " ha já tenha pago ", round(debt_paid), 
       "% de seu débito de extinção após ", 2026 - 1987, 
       " anos de isolamento.")

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
png("debt_paid_area_11.png", width = 14, height = 10,units = "cm", res = 900)
plot(debt_paid ~ log10(richness_2026$A),
     main = "Paid extinction debt over 39 years of island creation",
     xlab = "Island area (ha)",
     ylab = "Extinction debt paid (%)",
     xaxt = "n",
     las = 1,
     pch = 21, col = "black", bg = "lightgrey", cex = 1.5)+
  axis(1,
       at = log10(c(1, 10, 100, 1000)),
       labels = c(1, 10, 100, 1000))+
  abline(v = log10(c(1, 10, 100, 1000)),
         lty = "dashed", col = "lightgrey")
dev.off()


# Mean half-life between islands
cat("Mean half-life between islands:",
    mean(ilha_params$t_half_life), "\n")
max(ilha_params$t_half_life)
min(ilha_params$t_half_life)

t_half_ilhas <- ilha_params %>%
  select(ilhas, A, k0_ilha, k_ilha, t_half_life)
print(t_half_ilhas)

estimar_meia_vida <- function(A) {
  k <- k0_coef * A^alpha_coef
  t_half <- k * log(2)
  data.frame(A = A, k = k, t_half_life = t_half)}

estimar_meia_vida(100)


# Relaxation curves by island
t_seq <- seq(0, 100, length.out = 200)

curvas <- ilha_params %>%
  crossing(t = t_seq) %>%
  mutate(S_eq = c_coef * A^z_coef,
         S_pred = S_eq + (s_prevista - S_eq) * exp(-t / k_ilha))

biogeograp_model_plot <-
  ggplot() +
  geom_line(data = curvas,
            aes(x = t, y = S_pred, group = ilhas, colour = A),
            linewidth = 1) +
  geom_point(
    data = data_debt,
    aes(x = t, y = obs.richness, colour = A),
    size = 2, alpha = 0.6,
    position = position_jitter(width = 2, height = 0)) +
  
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  
  scale_colour_viridis_c(name = "Island area (ha)") +
  
  labs(x = "Time since isolation (years)",
       y = "Number of species remaining") +
  
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size = 10))
biogeograp_model_plot

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("11_biogeograp_model.png",
       plot = biogeograp_model_plot,
       width = 24,
       height = 12,
       units = "cm",
       dpi = 900)


####################################
##### Plot with reference areas ####

# Indicating the reference values for area
areas_ref <- c(1, 10, 100, 1000)
s_prevista_ref <- 10^(predict(m0, newdata = data.frame(a = areas_ref)))

# Predicted curves
curvas_ref <- expand.grid(A = areas_ref, t = t_seq) %>%
  mutate(
    S_eq   = c_coef * A^z_coef,
    k      = k0_coef * A^alpha_coef,
    S_pred = S_eq + (s_prevista_ref - S_eq) * exp(-t / k),
    A = factor(A, levels = areas_ref))

biogeograp_model_plot_ref <-
  ggplot() +
  geom_line(data = curvas,
            aes(x = t, y = S_pred, group = ilhas),
            colour = "grey85", linewidth = 0.4) +
  # geom_point(data = data_debt, aes(x = t, y = obs.richness),
  #           colour = "black", size = 2, alpha = 0.4) +
  geom_line(data = curvas_ref,
            aes(x = t, y = S_pred, group = A, colour = A),
            linewidth = 1) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  scale_colour_viridis_d(name = "Island area (ha)") +
  labs(x = "Time since isolation (years)",
       y = "Number of species") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size = 10))

biogeograp_model_plot_ref
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave("11_biogeograp_model_reference.png",
       plot = biogeograp_model_plot_ref,
       width = 24,
       height = 12,
       units = "cm",
       dpi = 900)
