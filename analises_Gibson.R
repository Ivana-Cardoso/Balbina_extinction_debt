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

gibson_island$t <- gibson_island$year - 1987
print(table(gibson_island$year, gibson_island$t))

gibson_SA <- gibson_island[,c(3,16)]
gibson_SA$area <- as.numeric(gibson_SA$area)
gibson_SA$richness <- as.numeric(gibson_SA$richness)
gibson_SA <- as.data.frame(gibson_SA)
SAR_MOD <- sar_power(gibson_SA)
SAR_MOD

fit <- nls(richness ~ sinf - (sinf - c*area^z) * exp(-k*t),
           data = gibson_island,
           start = list(sinf = 0.1, c = 1.1541884, z = 0.4214402, k = 0.1),
           control = nls.control(maxiter = 200))

print(summary(fit))
cat("\nCoeficientes:\n"); print(coef(fit))

summary(fit)


# R^2 (pseudo, como no artigo: 1 - SSresid/SStotal)
pred <- predict(fit)
ss_res <- sum((gibson_island$richness - pred)^2)
ss_tot <- sum((gibson_island$richness - mean(gibson_island$richness))^2)
r2 <- 1 - ss_res/ss_tot
cat("\nR^2:", r2, "\n")
