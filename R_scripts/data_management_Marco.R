
#### Marco 2010 data ####
rm(list = ls())
gc()
options(scipen = 999)
set.seed(13)

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco/Dados_originais")

# Import bird capture data
field_data <- read.csv("Marco_understory_birds_captures.csv", header = TRUE)

# Exclude recaptured individuals
field_data <- subset(
  field_data,
  recaptura == "N"
)

# Fix the species names
field_data$Espécie <- gsub(" ", "_", field_data$Espécie)
field_data$Espécie[field_data$Espécie == "Glaucidium_hardyii"] <- "Glaucidium_hardyi"
field_data$Espécie[field_data$Espécie == "Leptotila_rufaxila"] <- "Leptotila_rufaxilla"
field_data$Espécie[field_data$Espécie == "Willisornis_poecilonotus"] <- "Willisornis_poecilinotus"
field_data$Espécie[field_data$Espécie == "Xyphorhynchus_pardalotus"] <- "Xiphorhynchus_pardalotus"

# Create occurrence and unique ID column, sum the transects that are on the same island
colnames(field_data)[1] <- "Local"
colnames(field_data)[4] <- "Data"

field_data$Data <- as.Date(field_data$Data, format = "%d/%m/%Y")

# Some of Marco's islands were sampled by Anderson and me, so I will add a new 
# column with our islands names
field_data$Local[field_data$Local == "213 B"] <- "Cipoal"
field_data$Local[field_data$Local == "13"] <- "Piquia"
field_data$Local[field_data$Local == "18.4"] <- "Coata"
field_data$Local[field_data$Local == "475 C"] <- "Martelo"
field_data$Local[field_data$Local == "475 D"] <- "Martelo"
field_data$Local[field_data$Local == "1815 B"] <- "Gaviao_real"
field_data$Local[field_data$Local == "1815 C"] <- "Gaviao_real"
field_data$Local[field_data$Local == "16.5"] <- "Abusado"
field_data$Local[field_data$Local == "4.7"] <- "Panema"
field_data$Local[field_data$Local == "690 B"] <- "Beco_do_catitu"
field_data$Local[field_data$Local == "690 C"] <- "Beco_do_catitu"
field_data$Local[field_data$Local == "690 D"] <- "Beco_do_catitu"
field_data$Local[field_data$Local == "126"] <- "Pontal"
field_data$Local[field_data$Local == "85"] <- "Relogio"
field_data$Local[field_data$Local == "98.6"] <- "Sapupara"
field_data$Local[field_data$Local == "FCA LOR A"] <- "CF_Loreno"
field_data$Local[field_data$Local == "FCC WB A"] <- "CF_WABA"
field_data$Local[field_data$Local == "FCC WB B"] <- "CF_WABA"

field_data <- field_data %>%
  mutate(
    occurrence = 1,
    Ilha_Ano = paste(Local, Ano, sep = "_")
  )

# Construct a site-by-species matrix for island sites
comm <- field_data %>%
  group_by(Ilha_Ano, Local, Ano, Data, Espécie) %>%
  summarise(
    occurrence = sum(occurrence),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Espécie,
    values_from = occurrence,
    values_fill = 0
  )

insularization_year <- "01/10/1987"
insularization_year <- as.Date(insularization_year, format = "%d/%m/%Y")

comm$t <- time_length(
  interval(insularization_year, comm$Data),
  unit = "years"
)

comm$t <- round(comm$t, 2)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
#write.csv(comm, "comm_Marco.csv")



#### Bueno 2015 data ####
rm(list = ls())
gc()

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco/Dados_originais")

# Import bird capture data
field_data <- read.csv("Bueno_understory_birds_captures.csv", header = TRUE)

# Exclude recaptured individuals
field_data <- subset(
  field_data,
  status != "Recapture"
)

# Create occurrence and unique ID column, sum the transects that are on the same island
colnames(field_data)[15] <- "Local"
colnames(field_data)[11] <- "Data"
colnames(field_data)[13] <- "Ano"

field_data$Data <- as.Date(field_data$Data, format = "%d/%m/%Y")

field_data <- field_data %>%
  mutate(
    occurrence = 1,
    Ilha_Ano = paste(Local, Ano, sep = "_")
  )

# Construct a site-by-species matrix for island sites
comm <- field_data %>%
  group_by(Ilha_Ano, Local, Ano, Data, species) %>%
  summarise(
    occurrence = sum(occurrence),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = species,
    values_from = occurrence,
    values_fill = 0
  )

insularization_year <- "01/10/1987"
insularization_year <- as.Date(insularization_year, format = "%d/%m/%Y")

comm$t <- time_length(
  interval(insularization_year, comm$Data),
  unit = "years"
)

comm$t <- round(comm$t, 2)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
#write.csv(comm, "comm_Bueno.csv")


#### Amarante 2023 data ####
rm(list = ls())
gc()

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco/Dados_originais")

# Import bird capture data
field_data <- read.csv("Amarante_understory_birds_captures.csv", header = TRUE)

# Exclude recaptured individuals
field_data <- subset(
  field_data,
  status != "Recapture"
)

# Create occurrence and unique ID column, sum the transects that are on the same island
colnames(field_data)[15] <- "Local"
colnames(field_data)[11] <- "Data"
colnames(field_data)[13] <- "Ano"

field_data$Data <- as.Date(field_data$Data, format = "%d/%m/%Y")

field_data <- field_data %>%
  mutate(
    occurrence = 1,
    Ilha_Ano = paste(Local, Ano, sep = "_")
  )

# Construct a site-by-species matrix for island sites
comm <- field_data %>%
  group_by(Ilha_Ano, Local, Ano, Data, species) %>%
  summarise(
    occurrence = sum(occurrence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(species)) %>%
  pivot_wider(
    names_from = species,
    values_from = occurrence,
    values_fill = 0
  )
# note that it excluded the islands that were sampled, but no species were captured. I will add it back manually in excel

insularization_year <- "01/10/1987"
insularization_year <- as.Date(insularization_year, format = "%d/%m/%Y")

arrepiado_2024 <- "21/06/2024"
arrepiado_2024 <- as.Date(arrepiado_2024, format = "%d/%m/%Y")

time_length(
  interval(insularization_year, arrepiado_2024),
  unit = "years"
)

comm$t <- time_length(
  interval(insularization_year, comm$Data),
  unit = "years"
)

comm$t <- round(comm$t, 2)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
#write.csv(comm, "comm_Amarante.csv")
