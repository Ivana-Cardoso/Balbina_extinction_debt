#### Analises com os dados do Gibson ####

rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)

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
