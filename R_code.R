# Extinction debt in Balbina
#
# Ivana Cardoso
# ivanawaters@gmail.com
# Created on May 27, 2026
# Last modified on May 27, 2026

# I can't analyze the data separately by year (e.g., 2010, 2011, 2015, or 2016) because Marco did not resample any of the islands. He sampled each island only once, either in 2010 or in 2011.

library(dplyr)
library(iNEXT)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(sars)

# Preparing workspace
rm(list = ls())
gc()
options(scipen = 999)
set.seed(13)

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")

#### Preparing data ####
# Import bird capture data
Aurelio_Silva_comm <- read.csv("comm_Marco.csv", header = TRUE)
Bueno_comm <- read.csv("comm_Bueno.csv", header = TRUE)
Amarante_comm <- read.csv("comm_Amarante.csv", header = TRUE)

# Remove the continuous forest sites
Aurelio_Silva_comm <- Aurelio_Silva_comm[-c(11:12),]
Aurelio_Silva_comm <- bind_cols(
  Aurelio_Silva_comm[, 1:6],
  Aurelio_Silva_comm[, 7:105] %>% select(where(~ sum(.x, na.rm = TRUE) > 0))
)

Bueno_comm <- Bueno_comm[-c(13:22),]
Bueno_comm <- bind_cols(
  Bueno_comm[, 1:6],
  Bueno_comm[, 7:136] %>% select(where(~ sum(.x, na.rm = TRUE) > 0))
)

Amarante_comm <- Amarante_comm[-c(13:22),]
Amarante_comm <- bind_cols(
  Amarante_comm[, 1:6],
  Amarante_comm[, 7:120] %>% select(where(~ sum(.x, na.rm = TRUE) > 0))
)

dados <- bind_rows(Aurelio_Silva_comm, Bueno_comm, Amarante_comm)
dados[is.na(dados)] <- 0
colnames(dados)[1] <- "sites_6year"

dados$periodo_3 <- NA
dados$periodo_3[1:21] <- "2010"
dados$periodo_3[22:87] <- "2015"
dados$periodo_3[88:151] <- "2023"
dados$periodo_3 <- as.numeric(dados$periodo_3)

dados$t3 <- NA
dados$t3 <- dados$periodo_3 - 1987

dados$t5 <- NA
dados$t5 <- dados$year - 1987

dados$site_3year <- paste(dados$site, dados$periodo_3, sep = "_")

col_sp <- colnames(dados)[7:158]
dados_3years <- dados %>%
  group_by(site_3year) %>%
  summarise(
    area = first(area),
    t3 = first(t3),
    across(all_of(col_sp), sum, na.rm = TRUE)
  )

##############################################################
#### USING 6 PERIODS - 2010, 2011, 2015, 2016, 2023, 2014 ####
##############################################################

#### using iNEXT for richness (q=0) ####
build_abund_list <- function(df, sp_cols) {
  mat <- t(as.matrix(df[, sp_cols]))
  colnames(mat) <- df$site_year
  setNames(
    lapply(colnames(mat), function(nm) {
      v <- mat[, nm]
      v[v > 0 & !is.na(v)]      # remove especies ausentes naquela ilha
    }),
    colnames(mat)
  )
}

sp_marco    <- names(Aurelio_Silva_comm)[7:ncol(Aurelio_Silva_comm)]
sp_bueno    <- names(Bueno_comm)[7:ncol(Bueno_comm)]
sp_amarante <- names(Amarante_comm)[7:ncol(Amarante_comm)]

list_2010 <- build_abund_list(Aurelio_Silva_comm, sp_marco)
list_2015 <- build_abund_list(Bueno_comm, sp_bueno)
list_2023 <- build_abund_list(Amarante_comm, sp_amarante)

dados_5 <- c(list_2010, list_2015, list_2023)
dados_5 <- dados_5[sapply(dados_5, sum) > 0] # removes two islands (Arrepiado2024 and Formiga 2024 because they had 0 spp) because with them iNEXT doesnt work

x <- iNEXT(dados_5, q = 0, datatype = "abundance", endpoint = NULL)

# Returning the SC for each Assemblage after doubling the number of individuals.
SC <- aggregate(x$iNextEst$coverage_based$SC, 
          list(x$iNextEst$coverage_based$Assemblage), max)[, "x"]
hist(SC)
table(SC >= 0.8)

out_0.8 <- estimateD(dados_5, q = 0, datatype = "abundance",
                 base = "coverage", level = 0.8,
                 nboot = 1000,
                 conf = 0.95)
colnames(out_0.8)[1] <- "sites_6year"

out_0.9 <- estimateD(dados_5, q = 0, datatype = "abundance",
                 base = "coverage", level = 0.9,
                 nboot = 1000,
                 conf = 0.95)
colnames(out_0.9)[1] <- "sites_6year"

chao1 <- ChaoRichness(dados_5, datatype = "abundance")
chao1$sites_6year <- row.names(chao1)

dados <- dados %>%
  left_join(out_0.8 %>% select(sites_6year, qD), by = "sites_6year")
colnames(dados)[163]<- "est_rich_0.8"

dados <- dados %>%
  left_join(out_0.9 %>% select(sites_6year, qD), by = "sites_6year")
colnames(dados)[164] <- "est_rich_0.9"

dados <- dados %>%
  left_join(chao1 %>% select(sites_6year, Estimator), by = "sites_6year")
colnames(dados)[165] <- "chao1"

dados <- dados %>%
  left_join(chao1 %>% select(sites_6year, Observed), by = "sites_6year")
colnames(dados)[166] <- "obs.richness"

dados$obs.richness[is.na(dados$obs.richness)] <- 0

dados_6years <- dados

# setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt")
# write.csv(dados, "dados_completos_6periodos.csv")
rownames(dados_6years) <- dados_6years$sites_6year

##### Calculating the biogeographic model ####
# Following Gibson et al., 2013
# St = s_inf - (s_inf - c * a^z) * exp(-k * t)
#
# Parameters:
#   s_inf : asymptotic species richness as t -> infinity
#   c     : SAR constant  (S0 = c * a^z)
#   z     : SAR exponent
#   k     : relaxation rate (= I0 + E0)

########################################
#### Observed richness x 6 periodos ####
########################################
data_SA <- dados_6years[,c(6, 166)] #area and richness obs
data_SA$area <- as.numeric(data_SA$area)
data_SA$obs.richness <- as.numeric(data_SA$obs.richness)
data_SA <- as.data.frame(data_SA)
SAR_MOD <- sar_power(data_SA)
SAR_MOD # c = 3.6431114, z = 0.2382427

data_SA$t <- dados_6years$t5

fit <- nls(obs.richness ~ sinf - (sinf - c*area^z) * exp(-k*t),
                       data = data_SA,
                       start = list(sinf = 1, c = SAR_MOD$par[1], z = SAR_MOD$par[2], 
                                    k = 0.01),
                       control = nls.control(maxiter = 1000))

print(summary(fit))
cat("\nCoeficientes:\n"); print(coef(fit))


# R^2 (pseudo, como no artigo: 1 - SSresid/SStotal)
obs <- data_SA$obs.richness
pred <- predict(fit)
R2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
R2 #0.668

# Parâmetros 
s_inf <- coef(fit)[1]
cc <- coef(fit)[2]
zz <- coef(fit)[3]
k <- coef(fit)[4]

St <- function(tempo, area) {
  S0 <- cc * area^zz
  s_inf - (s_inf - S0) * exp(-k * tempo)
}

# derivada = taxa de extinção (painel B)
dSt <- function(tempo, area) {
  S0 <- cc * area^zz
  -k * (S0 - s_inf) * exp(-k * tempo)
}

# tempo até perder metade das espécies (painel C)
t_half <- function(area) {
  S0 <- cc * area^zz
  -(1/k) * log((s_inf - S0/2) / (s_inf - S0))
}

t_seq <- seq(0, 100, length.out = 300)

png("obs_rich_6t.png", width=1200, height=400, res=120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1)) 

areas_ilustrativas <- c(1, 5, 10, 25, 50, 100, 500, 1000)
cols <- c("cyan3","blue","darkgreen","red","black", "orange", "yellow", "pink")

plot(NULL, xlim=c(0,100), ylim=c(0,70),
     xlab="Time since isolation", ylab="Number of species remaining")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, St(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("topright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


plot(NULL, xlim=c(0,100), ylim=c(-1,0),
     xlab="Time since isolation", ylab="Rate of species extinction")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, dSt(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("bottomright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


areas_reais <- unique(data_SA$area)
areas_reais <- sort(areas_reais)

th <- t_half(areas_reais)
mean(th)
sd(th)

plot(areas_reais, th, type = "n",
     xlab = "Fragment area (ha)", ylab = expression(t[1/2]),
     xlim = c(0, 1850), ylim = c(17, 22))

# só plota/liga os pontos com t1/2 válido (a partir de ~0.8 ha)
ok <- !is.na(th)
points(areas_reais[ok], th[ok], pch = 1)
lines(areas_reais[ok], th[ok])

dev.off()
par(mfrow = c(1, 1))


################################################
#### Estimated richness SC 0.8 x 6 periodos ####
################################################
data_SA <- dados_6years[,c(6, 163)] 
data_SA$area <- as.numeric(data_SA$area)
data_SA$est_rich_0.8 <- as.numeric(data_SA$est_rich_0.8)
data_SA <- as.data.frame(data_SA)
data_SA <- na.omit(data_SA)
SAR_MOD <- sar_power(data_SA)
SAR_MOD # c = 4.6638544, z = 0.2350502

data_SA$t <- dados_6years$t5[-c(95,105)]

fit <- nls(est_rich_0.8 ~ sinf - (sinf - c*area^z) * exp(-k*t),
           data = data_SA,
           start = list(sinf = 1, c = 4.6638544, z = 0.2350502, 
                        k = 0.01),
           control = nls.control(maxiter = 1000))

print(summary(fit))
cat("\nCoeficientes:\n"); print(coef(fit))


# R^2 (pseudo, como no artigo: 1 - SSresid/SStotal)
obs <- data_SA$est_rich_0.8
pred <- predict(fit)
R2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
R2 #0.449

# Parâmetros 
s_inf <- coef(fit)[1]
cc <- coef(fit)[2]
zz <- coef(fit)[3]
k <- coef(fit)[4]

St <- function(tempo, area) {
  S0 <- cc * area^zz
  s_inf - (s_inf - S0) * exp(-k * tempo)
}

# derivada = taxa de extinção (painel B)
dSt <- function(tempo, area) {
  S0 <- cc * area^zz
  -k * (S0 - s_inf) * exp(-k * tempo)
}

# tempo até perder metade das espécies (painel C)
t_half <- function(area) {
  S0 <- cc * area^zz
  -(1/k) * log((s_inf - S0/2) / (s_inf - S0))
}

t_seq <- seq(0, 100, length.out = 300)

png("est_rich08_5t.png", width=1200, height=400, res=120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1)) 

areas_ilustrativas <- c(1, 5, 10, 25, 50, 100, 500, 1000)   # números redondos, só p/ ilustrar
cols <- c("cyan3","blue","darkgreen","red","black", "orange", "yellow", "pink")

plot(NULL, xlim=c(0,100), ylim=c(0,70),
     xlab="Time since isolation", ylab="Number of species remaining")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, St(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("topright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


plot(NULL, xlim=c(0,100), ylim=c(-1,0),
     xlab="Time since isolation", ylab="Rate of species extinction")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, dSt(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("bottomright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


areas_reais <- unique(data_SA$area)
areas_reais <- sort(areas_reais)

th <- t_half(areas_reais)
mean(th)
sd(th)

plot(areas_reais, th, type = "n",
     xlab = "Fragment area (ha)", ylab = expression(t[1/2]),
     xlim = c(0, 1850), ylim = c(20, 35))

# só plota/liga os pontos com t1/2 válido (a partir de ~0.8 ha)
ok <- !is.na(th)
points(areas_reais[ok], th[ok], pch = 1)
lines(areas_reais[ok], th[ok])

dev.off()
par(mfrow = c(1, 1))


############################
#### Chao1 x 6 periodos ####
############################
data_SA <- dados_6years[,c(6, 165)] 
data_SA$area <- as.numeric(data_SA$area)
data_SA$chao1 <- as.numeric(data_SA$chao1)
data_SA <- as.data.frame(data_SA)
data_SA <- na.omit(data_SA)
SAR_MOD <- sar_power(data_SA)
SAR_MOD # c = 6.9165520, z = 0.2461852

data_SA$t <- dados_6years$t5[-c(95,105)]

fit <- nls(chao1 ~ sinf - (sinf - c*area^z) * exp(-k*t),
           data = data_SA,
           start = list(sinf = 1, c = SAR_MOD$par[1], z = SAR_MOD$par[2], 
                        k = 0.01),
           control = nls.control(maxiter = 1000))

print(summary(fit))
cat("\nCoeficientes:\n"); print(coef(fit))


# R^2 (pseudo, como no artigo: 1 - SSresid/SStotal)
obs <- data_SA$chao1
pred <- predict(fit)
R2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
R2

# Parâmetros 
s_inf <- coef(fit)[1]
cc <- coef(fit)[2]
zz <- coef(fit)[3]
k <- coef(fit)[4]

St <- function(tempo, area) {
  S0 <- cc * area^zz
  s_inf - (s_inf - S0) * exp(-k * tempo)
}

# derivada = taxa de extinção (painel B)
dSt <- function(tempo, area) {
  S0 <- cc * area^zz
  -k * (S0 - s_inf) * exp(-k * tempo)
}

# tempo até perder metade das espécies (painel C)
t_half <- function(area) {
  S0 <- cc * area^zz
  -(1/k) * log((s_inf - S0/2) / (s_inf - S0))
}

t_seq <- seq(0, 100, length.out = 300)

png("chao1_5t.png", width=1200, height=400, res=120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1)) 

areas_ilustrativas <- c(1, 5, 10, 25, 50, 100, 500, 1000)   # números redondos, só p/ ilustrar
cols <- c("cyan3","blue","darkgreen","red","black", "orange", "yellow", "pink")

plot(NULL, xlim=c(0,100), ylim=c(0,100),
     xlab="Time since isolation", ylab="Number of species remaining")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, St(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("topright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda

plot(NULL, xlim=c(0,100), ylim=c(-1,0),
     xlab="Time since isolation", ylab="Rate of species extinction")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, dSt(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("bottomright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


areas_reais <- unique(data_SA$area)
areas_reais <- sort(areas_reais)

th <- t_half(areas_reais)
mean(th)
sd(th)

plot(areas_reais, th, type = "n",
     xlab = "Fragment area (ha)", ylab = expression(t[1/2]),
     xlim = c(0, 1850), ylim = c(14, 23))

ok <- !is.na(th)
points(areas_reais[ok], th[ok], pch = 1)
lines(areas_reais[ok], th[ok])

dev.off()
par(mfrow = c(1, 1))






























################################################################################
#### USING 3 PERIODS - 2010 (2010+2011), 2015 (2015+2016), 2023 (2023+2024) ####
################################################################################
Aurelio_Silva_comm$periodo_3 <- NA
Aurelio_Silva_comm$periodo_3 <- "2010"
Aurelio_Silva_comm$periodo_3 <- as.numeric(Aurelio_Silva_comm$periodo_3)

Aurelio_Silva_comm$t3 <- NA
Aurelio_Silva_comm$t3 <- Aurelio_Silva_comm$periodo_3 - 1987

Aurelio_Silva_comm$site_3year <- paste(Aurelio_Silva_comm$site, Aurelio_Silva_comm$periodo_3, sep = "_")

col_sp <- colnames(Aurelio_Silva_comm)[7:84]
Aurelio_Silva_3 <- Aurelio_Silva_comm %>%
  group_by(site_3year) %>%
  summarise(
    area = first(area),
    t3 = first(t3),
    across(all_of(col_sp), sum, na.rm = TRUE)
  )

Bueno_comm$periodo_3 <- NA
Bueno_comm$periodo_3 <- "2015"
Bueno_comm$periodo_3 <- as.numeric(Bueno_comm$periodo_3)

Bueno_comm$t3 <- NA
Bueno_comm$t3 <- Bueno_comm$periodo_3 - 1987

Bueno_comm$site_3year <- paste(Bueno_comm$site, Bueno_comm$periodo_3, sep = "_")

col_sp <- colnames(Bueno_comm)[7:115]
Bueno_comm_3 <- Bueno_comm %>%
  group_by(site_3year) %>%
  summarise(
    area = first(area),
    t3 = first(t3),
    across(all_of(col_sp), sum, na.rm = TRUE)
  )

Amarante_comm$periodo_3 <- NA
Amarante_comm$periodo_3 <- "2023"
Amarante_comm$periodo_3 <- as.numeric(Amarante_comm$periodo_3)

Amarante_comm$t3 <- NA
Amarante_comm$t3 <- Amarante_comm$periodo_3 - 1987

Amarante_comm$site_3year <- paste(Amarante_comm$site, Amarante_comm$periodo_3, sep = "_")

col_sp <- colnames(Amarante_comm)[7:95]
Amarante_comm_3 <- Amarante_comm %>%
  group_by(site_3year) %>%
  summarise(
    area = first(area),
    t3 = first(t3),
    across(all_of(col_sp), sum, na.rm = TRUE)
  )


#### using iNEXT for richness (q=0) ####
build_abund_list <- function(df, sp_cols) {
  mat <- t(as.matrix(df[, sp_cols]))
  colnames(mat) <- df$site_3year
  setNames(
    lapply(colnames(mat), function(nm) {
      v <- mat[, nm]
      v[v > 0 & !is.na(v)]      # remove especies ausentes naquela ilha
    }),
    colnames(mat)
  )
}
sp_marco    <- names(Aurelio_Silva_3)[7:ncol(Aurelio_Silva_3)]
sp_bueno    <- names(Bueno_comm_3)[7:ncol(Bueno_comm_3)]
sp_amarante <- names(Amarante_comm_3)[7:ncol(Amarante_comm_3)]

list_2010 <- build_abund_list(Aurelio_Silva_3, sp_marco)
list_2015 <- build_abund_list(Bueno_comm_3, sp_bueno)
list_2023 <- build_abund_list(Amarante_comm_3, sp_amarante)

dados_3 <- c(list_2010, list_2015, list_2023)
dados_3 <- dados_3[sapply(dados_3, sum) > 0]

x <- iNEXT(dados_3, q = 0, datatype = "abundance", endpoint = NULL)

# Returning the SC for each Assemblage after doubling the number of individuals.
SC <- aggregate(x$iNextEst$coverage_based$SC, 
                list(x$iNextEst$coverage_based$Assemblage), max)[, "x"]
hist(SC)
table(SC >= 0.8)

out_0.8 <- estimateD(dados_3, q = 0, datatype = "abundance",
                     base = "coverage", level = 0.8,
                     nboot = 1000,
                     conf = 0.95)
colnames(out_0.8)[1] <- "sites_3year"

out_0.9 <- estimateD(dados_3, q = 0, datatype = "abundance",
                     base = "coverage", level = 0.9,
                     nboot = 1000,
                     conf = 0.95)
colnames(out_0.9)[1] <- "sites_3year"

chao1 <- ChaoRichness(dados_3, datatype = "abundance")
chao1$sites_3year <- row.names(chao1)

colnames(dados_3years)[1] <- "sites_3year"

dados_3years <- dados_3years %>%
  left_join(out_0.8 %>% select(sites_3year, qD), by = "sites_3year")
colnames(dados_3years)[156]<- "est_rich_0.8"

dados_3years <- dados_3years %>%
  left_join(out_0.9 %>% select(sites_3year, qD), by = "sites_3year")
colnames(dados_3years)[157] <- "est_rich_0.9"

dados_3years <- dados_3years %>%
  left_join(chao1 %>% select(sites_3year, Estimator), by = "sites_3year")
colnames(dados_3years)[158] <- "chao1"

dados_3years <- dados_3years %>%
  left_join(chao1 %>% select(sites_3year, Observed), by = "sites_3year")
colnames(dados_3years)[159] <- "obs.richness"

dados_3years$obs.richness[is.na(dados_3years$obs.richness)] <- 0

# setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt")
# write.csv(dados_3years, "dados_completos_3periodos.csv")
rownames(dados_3years) <- dados_3years$sites_3year

##### Calculating the biogeographic model ####
# Following Gibson et al., 2013
# St = s_inf - (s_inf - c * a^z) * exp(-k * t)
#
# Parameters:
#   s_inf : asymptotic species richness as t -> infinity
#   c     : SAR constant  (S0 = c * a^z)
#   z     : SAR exponent
#   k     : relaxation rate (= I0 + E0)

########################################
#### Observed richness x 3 periodos ####
########################################
data_SA <- dados_3years[,c(2, 159)] #area and richness obs
data_SA$area <- as.numeric(data_SA$area)
data_SA$obs.richness <- as.numeric(data_SA$obs.richness)
data_SA <- as.data.frame(data_SA)
SAR_MOD <- sar_power(data_SA)
SAR_MOD # c = 4.3018665, z = 0.2545439

data_SA$t <- dados_3years$t3

fit <- nls(obs.richness ~ sinf - (sinf - c*area^z) * exp(-k*t),
           data = data_SA,
           start = list(sinf = 1, c = SAR_MOD$par[1], z = SAR_MOD$par[2], 
                        k = 0.01),
           control = nls.control(maxiter = 1000))

print(summary(fit))
cat("\nCoeficientes:\n"); print(coef(fit))


# R^2 (pseudo, como no artigo: 1 - SSresid/SStotal)
obs <- data_SA$obs.richness
pred <- predict(fit)
R2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
R2

# Parâmetros 
s_inf <- coef(fit)[1]
cc <- coef(fit)[2]
zz <- coef(fit)[3]
k <- coef(fit)[4]

St <- function(tempo, area) {
  S0 <- cc * area^zz
  s_inf - (s_inf - S0) * exp(-k * tempo)
}

# derivada = taxa de extinção (painel B)
dSt <- function(tempo, area) {
  S0 <- cc * area^zz
  -k * (S0 - s_inf) * exp(-k * tempo)
}

# tempo até perder metade das espécies (painel C)
t_half <- function(area) {
  S0 <- cc * area^zz
  -(1/k) * log((s_inf - S0/2) / (s_inf - S0))
}

t_seq <- seq(0, 100, length.out = 300)

png("obs_rich_3t.png", width=1200, height=400, res=120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1)) 

areas_ilustrativas <- c(1, 5, 10, 25, 50, 100, 500, 1000)   # números redondos, só p/ ilustrar
cols <- c("cyan3","blue","darkgreen","red","black", "orange", "yellow", "pink")

plot(NULL, xlim=c(0,100), ylim=c(0,35),
     xlab="Time since isolation", ylab="Number of species remaining")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, St(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("topright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda

plot(NULL, xlim=c(0,100), ylim=c(-0.4,0),
     xlab="Time since isolation", ylab="Rate of species extinction")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, dSt(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("bottomright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


areas_reais <- unique(data_SA$area)
areas_reais <- sort(areas_reais)

th <- t_half(areas_reais)
mean(na.omit(th))
sd(na.omit(th))


plot(areas_reais, th, type = "n",
     xlab = "Fragment area (ha)", ylab = expression(t[1/2]),
     xlim = c(0, 1850), ylim = c(70, 300))

ok <- !is.na(th)
points(areas_reais[ok], th[ok], pch = 1)
lines(areas_reais[ok], th[ok])

dev.off()
par(mfrow = c(1, 1))



################################################
#### Estimated richness SC 0.8 x 3 periodos ####
################################################
data_SA <- dados_3years[,c(2, 156)] #area and est rich
data_SA$area <- as.numeric(data_SA$area)
data_SA$est_rich_0.8 <- as.numeric(data_SA$est_rich_0.8)
data_SA <- as.data.frame(data_SA)
data_SA <- na.omit(data_SA)
SAR_MOD <- sar_power(data_SA)
SAR_MOD # c = 6.2861128, z = 0.1903162

data_SA$t <- dados_3years$t3

fit <- nls(est_rich_0.8 ~ sinf - (sinf - c*area^z) * exp(-k*t),
           data = data_SA,
           start = list(sinf = 1, c = SAR_MOD$par[1], z = SAR_MOD$par[2], 
                        k = 0.01),
           control = nls.control(maxiter = 1000))

print(summary(fit))
cat("\nCoeficientes:\n"); print(coef(fit))


# R^2 (pseudo, como no artigo: 1 - SSresid/SStotal)
obs <- data_SA$est_rich_0.8
pred <- predict(fit)
R2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
R2

# Parâmetros 
s_inf <- coef(fit)[1]
cc <- coef(fit)[2]
zz <- coef(fit)[3]
k <- coef(fit)[4]

St <- function(tempo, area) {
  S0 <- cc * area^zz
  s_inf - (s_inf - S0) * exp(-k * tempo)
}

# derivada = taxa de extinção (painel B)
dSt <- function(tempo, area) {
  S0 <- cc * area^zz
  -k * (S0 - s_inf) * exp(-k * tempo)
}

# tempo até perder metade das espécies (painel C)
t_half <- function(area) {
  S0 <- cc * area^zz
  -(1/k) * log((s_inf - S0/2) / (s_inf - S0))
}

t_seq <- seq(0, 100, length.out = 300)

png("est_rich08_3t.png", width=1200, height=400, res=120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1)) 

areas_ilustrativas <- c(1, 5, 10, 25, 50, 100, 500, 1000)   # números redondos, só p/ ilustrar
cols <- c("cyan3","blue","darkgreen","red","black", "orange", "yellow", "pink")

plot(NULL, xlim=c(0,100), ylim=c(0,70),
     xlab="Time since isolation", ylab="Number of species remaining")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, St(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("topright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


plot(NULL, xlim=c(0,100), ylim=c(-1,0),
     xlab="Time since isolation", ylab="Rate of species extinction")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, dSt(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("bottomright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


areas_reais <- unique(data_SA$area)
areas_reais <- sort(areas_reais)

th <- t_half(areas_reais)
mean(na.omit(th))
sd(na.omit(th))

plot(areas_reais, th, type = "n",
     xlab = "Fragment area (ha)", ylab = expression(t[1/2]),
     xlim = c(0, 1850), ylim = c(20, 72))

ok <- !is.na(th)
points(areas_reais[ok], th[ok], pch = 1)
lines(areas_reais[ok], th[ok])

dev.off()
par(mfrow = c(1, 1))


#############################
#### Chao 1 x 3 periodos ####
#############################
data_SA <- dados_3years[,c(2, 158)] #area and est rich
data_SA$area <- as.numeric(data_SA$area)
data_SA$chao1 <- as.numeric(data_SA$chao1)
data_SA <- as.data.frame(data_SA)
data_SA <- na.omit(data_SA)
SAR_MOD <- sar_power(data_SA)
SAR_MOD # c = 10.028509, z = .210072

data_SA$t <- dados_3years$t3

fit <- nls(chao1 ~ sinf - (sinf - c*area^z) * exp(-k*t),
           data = data_SA,
           start = list(sinf = 1, c = SAR_MOD$par[1], z = SAR_MOD$par[2], 
                        k = 0.01),
           control = nls.control(maxiter = 1000))

print(summary(fit))
cat("\nCoeficientes:\n"); print(coef(fit))


obs <- data_SA$chao1
pred <- predict(fit)
R2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
R2

# Parâmetros 
s_inf <- coef(fit)[1]
cc <- coef(fit)[2]
zz <- coef(fit)[3]
k <- coef(fit)[4]

St <- function(tempo, area) {
  S0 <- cc * area^zz
  s_inf - (s_inf - S0) * exp(-k * tempo)
}

# derivada = taxa de extinção (painel B)
dSt <- function(tempo, area) {
  S0 <- cc * area^zz
  -k * (S0 - s_inf) * exp(-k * tempo)
}

# tempo até perder metade das espécies (painel C)
t_half <- function(area) {
  S0 <- cc * area^zz
  -(1/k) * log((s_inf - S0/2) / (s_inf - S0))
}

t_seq <- seq(0, 100, length.out = 300)

png("est_rich08_3t.png", width=1200, height=400, res=120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1)) 

areas_ilustrativas <- c(1, 5, 10, 25, 50, 100, 500, 1000)   # números redondos, só p/ ilustrar
cols <- c("cyan3","blue","darkgreen","red","black", "orange", "yellow", "pink")

plot(NULL, xlim=c(0,100), ylim=c(5,80),
     xlab="Time since isolation", ylab="Number of species remaining")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, St(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("topright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   


plot(NULL, xlim=c(0,100), ylim=c(-1,0),
     xlab="Time since isolation", ylab="Rate of species extinction")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, dSt(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}
legend("bottomright",
       legend = paste(areas_ilustrativas, "ha"),
       col    = cols,
       lwd    = 2,
       cex    = 0.8,
       bty    = "n")   # bty="n" tira a caixa ao redor da legenda


areas_reais <- unique(data_SA$area)
areas_reais <- sort(areas_reais)

th <- t_half(areas_reais)
mean(na.omit((th)))
sd(na.omit(th))

plot(areas_reais, th, type = "n",
     xlab = "Fragment area (ha)", ylab = expression(t[1/2]),
     xlim = c(0, 1850), ylim = c(40, 160))

ok <- !is.na(th)
points(areas_reais[ok], th[ok], pch = 1)
lines(areas_reais[ok], th[ok])

dev.off()
par(mfrow = c(1, 1))



## boxplots 
#1. 3 períodos e Sobs (o gráfico do projeto da tese)
#2. 5 períodos e Sobs
#3. 3 períodos e Sraref
#4. 5 períodos e Sraref


library(ggplot2)

dados_6years$year <- as.factor(dados_6years$year)

obs_rich_5 <-
  ggplot(mapping = aes(x = year, y = obs.richness),
         data = dados_6years) +
  labs(x = "Ano de amostragem", y = "Riqueza observada",
       title = "Riqueza observada - seis períodos") +
  geom_boxplot(outliers = FALSE) +
  geom_point(size = 3, alpha = 0.3) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "bottom")
obs_rich_5

est_rich_5 <-
  ggplot(mapping = aes(x = year, y = est_rich_0.8),
         data = dados_ano) +
  labs(x = "Ano de amostragem", y = "Riqueza estimada SC 0.8",
       title = "Riqueza estimada - cinco períodos") +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitter(), size = 3, alpha = 0.3) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "bottom")
est_rich_5

obs_rich_3 <-
  ggplot(mapping = aes(x = periodo_3, y = obs_rich),
         data = dados_ano) +
  labs(x = "Ano de amostragem", y = "Riqueza observada",
       title = "Riqueza observada - três períodos") +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(
    labels = c(
      "2010" = "2010-2011",
      "2015" = "2015-2016",
      "2023" = "2023-2024"
    )
  )+
  geom_point(position = position_jitter(), size = 3, alpha = 0.3) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "bottom")
obs_rich_3

est_rich_3 <-
  ggplot(mapping = aes(x = periodo_3, y = est_rich_0.8),
         data = dados_ano) +
  labs(x = "Ano de amostragem", y = "Riqueza estimada SC 0.8",
       title = "Riqueza estimada - três períodos") +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(
    labels = c(
      "2010" = "2010-2011",
      "2015" = "2015-2016",
      "2023" = "2023-2024"
    )
  )+
  geom_point(position = position_jitter(), size = 3, alpha = 0.3) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.25),
        plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
        legend.position = "bottom")
est_rich_3

library(patchwork)

plot1 <- (obs_rich_5 | est_rich_5) /
  (obs_rich_3 | est_rich_3)

plot1

ggsave("boxplots_balbina.png",
  plot = plot1,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 900
)

