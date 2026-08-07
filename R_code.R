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

# Preparing workspace
rm(list = ls())
gc()
options(scipen = 999)
set.seed(13)

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Data")

#### Preparing data ####
# Import bird capture data
Aurelio_Silva_comm <- read.csv("comm_Marco.csv", header = TRUE)
colnames(Aurelio_Silva_comm)[2] <- "site"

Bueno_comm <- read.csv("comm_Anderson.csv", header = TRUE)

Amarante_comm <- read.csv("comm_Amarante.csv", header = TRUE)
colnames(Amarante_comm)[1] <- "site"
Amarante_comm$site[11] <- "CF_WABA"
Amarante_comm$site[18] <- "Gaviao-real"

# Keep only the 11 islands (+2 continuous forest) in common between all of us
Aurelio_Silva_comm <- Aurelio_Silva_comm[c(2,3,4,6,7,9,14,15,17,18,19,20,21),-1]
Aurelio_Silva_comm <- bind_cols(
  Aurelio_Silva_comm[, 1:2],
  Aurelio_Silva_comm[, 3:97] %>% select(where(~ sum(.x, na.rm = TRUE) > 0))
)
rownames(Aurelio_Silva_comm) <- Aurelio_Silva_comm$site

Bueno_comm <- Bueno_comm[Bueno_comm$site %in% Aurelio_Silva_comm$site,]
Bueno_comm <- bind_cols(
  Bueno_comm[, 1],
  Bueno_comm[, 2:131] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)),
  Bueno_comm[,132:137]
)
colnames(Bueno_comm)[1] <- "site"
rownames(Bueno_comm) <- Bueno_comm$site

Amarante_comm <- Amarante_comm[Amarante_comm$site %in% Aurelio_Silva_comm$site,]
Amarante_comm <- bind_cols(
  Amarante_comm[, 1],
  Amarante_comm[, 4:117] %>% select(where(~ sum(.x, na.rm = TRUE) > 0)),
  Amarante_comm[,118]
)
colnames(Amarante_comm)[1] <- "site"
rownames(Amarante_comm) <- Amarante_comm$site
colnames(Amarante_comm)[76] <- "area"


#### Rarefaction ####
# Because there is a difference in sampling effort between Marco's study and ours, I need to use rarefaction curves to standardize all samples to the same number of individuals.
set.seed(13)

# using iNEXT
# Creating a function to convert each matrix into a named list (island_period -> abundance vector)
matrix_to_list <- function(df, periodo) {
  mat <- t(as.matrix(df))          
  colnames(mat) <- rownames(df)    
  lst <- as.list(as.data.frame(mat))
  names(lst) <- paste(names(lst), periodo, sep = "_")
  lst
}

list_2010 <- matrix_to_list(Aurelio_Silva_comm[,c(3:93)], "2010")
list_2015 <- matrix_to_list(Bueno_comm[,c(2:84)], "2015")
list_2023 <- matrix_to_list(Amarante_comm[,c(2:75)], "2023")

dados <- c(list_2010, list_2015, list_2023)

# Run iNEXT for richness (q=0)
out <- iNEXT(dados, q = 0, datatype = "abundance")

# Extract the size-based estimates data frame and split island/period back into separate columns
df_plot <- out$iNextEst$size_based %>%
  mutate(
    Periodo = str_extract(Assemblage, "\\d{4}$"),        # 4 last numbers
    Ilha    = str_remove(Assemblage, "_\\d{4}$")          # remove "_" + 4 last numbers
  )

df_plot_rarefaction <- df_plot %>%
  filter(Method %in% c("Rarefaction", "Observed"))

ggplot(df_plot_rarefaction, aes(x = m, y = qD, color = Periodo, group = Assemblage)) +
  geom_line(
    data = df_plot_rarefaction %>% filter(Method == "Rarefaction")
  ) +
  geom_point(
    data = df_plot_rarefaction %>% filter(Method == "Observed"),
    size = 2
  ) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Periodo), alpha = 0.15, color = NA) +
  facet_wrap(~ Ilha) +
  labs(x = "Nº de indivíduos", y = "Riqueza estimada (q=0)") +
  theme_minimal()


## Minimum sampling effort (n ind) shared across the three years for each island
ilhas <- unique(df_plot$Ilha)

resultados <- list()

for (ilha in ilhas) {
  # select the three datasets (time periods) for this island from the 'dados' list
  nomes_ilha <- names(dados)[str_remove(names(dados), "_\\d{4}$") == ilha]
  subset_dados <- dados[nomes_ilha]
  
  # smallest observed sample size (N) across the time periods for this island
  menor_N_ilha <- min(sapply(subset_dados, sum))
  
  # standardized estimateD at the minimum sample size (N)
  est <- estimateD(subset_dados, datatype = "abundance",
                   base = "size", level = menor_N_ilha, q = 0)
  
  resultados[[ilha]] <- est
}

tabela_comparacao <- bind_rows(resultados, .id = "IlhaGroup") %>%
  mutate(
    Periodo = str_extract(Assemblage, "\\d{4}$"),
    Ilha    = str_remove(Assemblage, "_\\d{4}$")
  )

tabela_comparacao

tabela_simples <- tabela_comparacao %>%
  select(site = Ilha, richness = qD, year = Periodo)

tabela_simples <- tabela_simples %>%
  arrange(year, site)

tabela_simples <- tabela_simples %>%
  mutate(year = case_when(
    year == "2010" ~ "2010-2011",
    year == "2015" ~ "2015-2016",
    year == "2023" ~ "2023-2024",
    TRUE ~ year  # mantém qualquer outro valor inesperado sem alterar
  ))

tabela_simples$t <- NA
tabela_simples$t[1:13] <- 2010.5 - 1987
tabela_simples$t[14:26] <- 2015.5 - 1987
tabela_simples$t[27:39] <- 2023.5 - 1987

tabela_simples$area <- NA

Aurelio_Silva_comm <- Aurelio_Silva_comm %>%
  arrange(Aurelio_Silva_comm$site)
Aurelio_Silva_comm$site == tabela_simples$site[1:13]
tabela_simples$area[1:13] <- Aurelio_Silva_comm$area

Bueno_comm <- Bueno_comm %>%
  arrange(Bueno_comm$site)
Bueno_comm$site == tabela_simples$site[14:26]
tabela_simples$area[14:26] <- Bueno_comm$area

Amarante_comm <- Amarante_comm %>%
  arrange(Amarante_comm$site)
Amarante_comm$site == tabela_simples$site[27:39]
tabela_simples$area[27:39] <- Amarante_comm$area

tabela_simples <- tabela_simples %>%
  filter(!site %in% c("CF_Loreno", "CF_WABA"))

tabela_simples$area <- as.numeric(tabela_simples$area)

plot = 
  ggplot(tabela_simples, 
         aes(x=area, y=richness,
                        colour=year)) +
  scale_x_log10() +
  
  geom_smooth(method = "lm", se = FALSE, aes(fill = year), alpha = 0.3, linewidth = 1)+
  
  geom_point(size=4, shape = 19, alpha = 0.6)+
  
  labs(y = "Rarefied number of species",
       x = "Island area (ha) log10",
       color = "Year",
       fill = "Year")+
  
  theme_bw(base_size = 10) +
  
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "#252525"),
        axis.title = element_text(colour = "#252525", face = "bold", size = 14),
        axis.title.x = element_text(vjust = 0, size = 14),
        axis.title.y = element_text(vjust = 2, size = 14),
        axis.text = element_text(colour = "#252525"),
        axis.ticks = element_line(colour = "#252525", size = 0.25))
plot  

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Balbina_cap.3/Figures")
ggsave(
  plot = plot,
  filename = "rarefied_sp_area.png", dpi = 600,
  width = 12, height = 10, units = 'cm')
