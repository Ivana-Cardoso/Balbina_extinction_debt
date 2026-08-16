#### Analises com os dados do Gibson ####

rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)
library(sars)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3")
gibson_transec <- read.csv("gibson_data_.csv")

# criando as categorias de tamanho das ilhas
gibson_transec <- gibson_transec %>%
  mutate(
    size = case_when(
      area >= 0.3 & area <= 4.7 ~ "Small",
      area >= 10.1 & area <= 56.3 ~ "Large",
      TRUE ~ NA_character_
    )
  )


dados_periodo <- gibson_transec %>%
  mutate(periodo = case_when(
    year %in% c(1992, 1993, 1994) ~ "1992-1994",
    year %in% c(2012, 2013) ~ "2012-2013",
    TRUE ~ NA_character_
  )) 

small_islands <- subset(dados_periodo, size == "Small")
large_islands <- subset(dados_periodo, size == "Large")

# Large islands
LI <-
ggplot(large_islands, aes(x = periodo, y = richness)) +
  labs(y = "Richness",
       x = "Years after isolation",
       title = "Large islands")+
  scale_y_continuous(limits = c(0, 14),
                     breaks = seq(0, 14, 2)) +
  geom_hline(yintercept = 13, linetype = "dashed") +
  geom_boxplot() +
  geom_point(size = 2) +
  theme_classic()
LI

#Small islands
SI <-
ggplot(small_islands, aes(x = periodo, y = richness)) +
  labs(y = "Richness",
       x = "Years after isolation",
       title = "Small islands") +
  scale_y_continuous(limits = c(0, 14),
                     breaks = seq(0, 14, 2)) +
  geom_hline(yintercept = 13, linetype = "dashed") +
  geom_boxplot() +
  geom_point(size = 2) +
  theme_classic()
SI

ggsave("Gibson_Fig2.png",
       LI + SI,
       width = 10,
       height = 5,
       dpi = 300)


###############################
##### Modelo Biogeografico ####
# Material suplementar (Eq. S5):
# St = s_inf - (s_inf - c * a^z) * exp(-k * t)
#
# Parameters:
#   s_inf : asymptotic species richness as t -> infinity
#   c     : SAR constant  (S0 = c * a^z)
#   z     : SAR exponent
#   k     : relaxation rate (= I0 + E0)

# Derivando os valores de c e z
# os valores obtidos devem ser:
# S_inf = 0.751
# c = 2.223
# z = 0.482
# k = 0.0732

gibson_data <- read.csv("Table_S1_small_mammal_abundance.csv")

gibson_island <- gibson_data %>%
  group_by(year, island, area) %>%
  summarise(
    across(
      Callosciurus.caniceps:Tupaia.glis,
      sum,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    richness = rowSums(across(Callosciurus.caniceps:Tupaia.glis) > 0)
  )

abund <- t(as.matrix(gibson_island[, 4:15]))
out <- iNEXT(abund)



gibson_island$t <- gibson_island$year - 1987
print(table(gibson_island$year, gibson_island$t))

gibson_island <- gibson_island[!gibson_island$island %in% c("X1", "X2", "X3", "X4"), ]


gibson_SA <- gibson_island[,c(3,16)]
gibson_SA$area <- as.numeric(gibson_SA$area)
gibson_SA$richness <- as.numeric(gibson_SA$richness)
gibson_SA <- as.data.frame(gibson_SA)
SAR_MOD <- sar_power(gibson_SA)
SAR_MOD

#gibson_SA$t <- gibson_island$t
#gibson_SA$island <- gibson_island$island
#gibson_SA$year <- gibson_island$year

fit <- nls(richness ~ sinf - (sinf - c*area^z) * exp(-k*t),
           data = gibson_island,
           start = list(sinf = 0.1, c = 1.3484891, z = 0.4371562, k = 0.1),
           control = nls.control(maxiter = 200))

#fit2 <- nls(richness ~ sinf - (sinf - c*area^z) * exp(-k*t),
#           data = gibson_SA,
#           start = list(sinf = 0.1, c = 1.1541884, z = 0.4214402, k = 0.1),
#           control = nls.control(maxiter = 200))

print(summary(fit))

cat("\nCoeficientes:\n"); print(coef(fit))
summary(fit)


# R^2 (pseudo, como no artigo: 1 - SSresid/SStotal)
pred <- predict(fit)
ss_res <- sum((gibson_island$richness - pred)^2)
ss_tot <- sum((gibson_island$richness - mean(gibson_island$richness))^2)
r2 <- 1 - ss_res/ss_tot
cat("\nR^2:", r2, "\n")



# Parâmetros do artigo
s_inf <- 0.751
cc    <- 2.223
zz    <- 0.482
k     <- 0.0732

# Parâmetros alcançados pela função
(s_inf <- coef(fit)[1])
(cc    <- coef(fit)[2])
(zz    <- coef(fit)[3])
(k     <- coef(fit)[4])

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

t_seq <- seq(0, 80, length.out = 300)

png("figura_gibson.png", width=1200, height=400, res=120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1)) 
## ---- Painel A: número de espécies remanescentes ----
areas_ilustrativas <- c(1, 5, 10, 25, 50)   # números redondos, só p/ ilustrar
cols <- c("cyan3","blue","darkgreen","red","black")

plot(NULL, xlim=c(0,80), ylim=c(0,15),
     xlab="Time since isolation", ylab="Number of species remaining")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, St(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}

## ---- Painel B: taxa de extinção ----
plot(NULL, xlim=c(0,80), ylim=c(-1,0),
     xlab="Time since isolation", ylab="Rate of species extinction")
for (i in seq_along(areas_ilustrativas)) {
  lines(t_seq, dSt(t_seq, areas_ilustrativas[i]), col=cols[i], lwd=2)
}

## ---- Painel C: t1/2 vs área real das ilhas ----
# Parâmetros do artigo

# áreas reais das 16 ilhas do estudo
areas_reais <- c(0.3, 0.4, 0.8, 1.0, 1.1, 1.4, 1.7, 1.9,
                 4.7, 10.1, 10.4, 12.1, 21.2, 23.5, 24.4, 56.3)

th <- t_half(areas_reais)

## Painel C
plot(areas_reais, th, type = "n",
     xlab = "Fragment area (ha)", ylab = expression(t[1/2]),
     xlim = c(0, 60), ylim = c(10, 22))

# só plota/liga os pontos com t1/2 válido (a partir de ~0.8 ha)
ok <- !is.na(th)
points(areas_reais[ok], th[ok], pch = 1)
lines(areas_reais[ok], th[ok])

dev.off()
par(mfrow = c(1, 1))