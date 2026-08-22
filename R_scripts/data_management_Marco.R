
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

# Fix the species names - CBRO 2021
field_data$Espécie <- gsub(" ", "_", field_data$Espécie)
field_data$Espécie[field_data$Espécie == "Glaucidium_hardyii"] <- "Glaucidium_hardyi"
field_data$Espécie[field_data$Espécie == "Leptotila_rufaxila"] <- "Leptotila_rufaxilla"
field_data$Espécie[field_data$Espécie == "Willisornis_poecilonotus"] <- "Willisornis_poecilinotus"
field_data$Espécie[field_data$Espécie == "Xyphorhynchus_pardalotus"] <- "Xiphorhynchus_pardalotus"
field_data$Espécie[field_data$Espécie == "Dixiphia_pipra"] <- "Pseudopipra_pipra"
field_data$Espécie[field_data$Espécie == "Platyrinchus_platyrinchos"] <- "Platyrinchus_platyrhynchos"
field_data$Espécie[field_data$Espécie == "Lanio_surinamus"] <- "Maschalethraupis_surinamus"
field_data$Espécie[field_data$Espécie == "Myrmotherula_guttata"] <- "Isleria_guttata"
field_data$Espécie[field_data$Espécie == "Tangara_palmarum"] <- "Thraupis_palmarum"
field_data$Espécie[field_data$Espécie == "Myrmeciza_ferruginea"] <- "Myrmoderus_ferrugineus"
field_data$Espécie[field_data$Espécie == "Tachyphonus_surinamus"] <- "Maschalethraupis_surinamus"
field_data$Espécie[field_data$Espécie == "Schistocichla_leucostigma"] <- "Myrmelastes_leucostigma"
field_data$Espécie[field_data$Espécie == "Cyanocompsa_cyanoides"] <- "Cyanoloxia_rothschildii"
field_data$Espécie[field_data$Espécie == "Deconychura_stictolaema"] <- "Certhiasomus_stictolaemus"
field_data$Espécie[field_data$Espécie == "Corapipo_guturallis"] <- "Corapipo_gutturalis"
field_data$Espécie[field_data$Espécie == "Myrmotherula_gutturalis"] <- "Epinecrophylla_gutturalis"
field_data$Espécie[field_data$Espécie == "Cercomacra_laeta"] <- "Cercomacroides_laeta"
field_data$Espécie[field_data$Espécie == "Cercomacra_tyrannina"] <- "Cercomacroides_tyrannina"
field_data$Espécie[field_data$Espécie == "Philydor_erythocercum"] <- "Philydor_erythrocercum"

# Create occurrence and unique ID column, sum the transects that are on the same island
colnames(field_data)[1] <- "Local"
colnames(field_data)[4] <- "Data"

field_data$Data <- as.Date(field_data$Data, format = "%d/%m/%Y")

# Some of Marco's islands were sampled by Anderson and me, so I will add a new 
# column with our islands names
field_data$Local[field_data$Local == "213 B"] <- "Cipoal"
field_data$Local[field_data$Local == "13"] <- "Piquia"
field_data$Local[field_data$Local == "18.4"] <- "Coata"
field_data$Local[field_data$Local == "475 C"] <- "Martelo_C"
field_data$Local[field_data$Local == "475 D"] <- "Martelo_D"
field_data$Local[field_data$Local == "1815 B"] <- "Gaviao_real_B"
field_data$Local[field_data$Local == "1815 C"] <- "Gaviao_real_C"
field_data$Local[field_data$Local == "16.5"] <- "Abusado"
field_data$Local[field_data$Local == "4.7"] <- "Panema"
field_data$Local[field_data$Local == "690 B"] <- "Beco_do_catitu_B"
field_data$Local[field_data$Local == "690 C"] <- "Beco_do_catitu_C"
field_data$Local[field_data$Local == "690 D"] <- "Beco_do_catitu_D"
field_data$Local[field_data$Local == "126"] <- "Pontal"
field_data$Local[field_data$Local == "85"] <- "Relogio"
field_data$Local[field_data$Local == "98.6"] <- "Sapupara"
field_data$Local[field_data$Local == "FCA LOR A"] <- "CF_Loreno"
field_data$Local[field_data$Local == "FCC WB A"] <- "CF_WABA_A"
field_data$Local[field_data$Local == "FCC WB B"] <- "CF_WABA_B"

field_data <- field_data %>%
  mutate(
    occurrence = 1,
    Ilha_Ano = paste(Local, Ano, sep = "_"))

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
comm <- comm %>%
  mutate(ilha = sub("_[A-Za-z]$", "", Local))

comm$t_period <- 2010-1987

comm <- comm[,c(1,2,103,3,4,102,104,5:101)]
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
write.csv(comm, "comm_Marco.csv")



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

field_data$species[field_data$species == "Amazilia_versicolor"] <- "Chrysuronia_versicolor"
field_data$species[field_data$species == "Dixiphia_pipra"] <- "Pseudopipra_pipra"
field_data$species[field_data$species == "Lanio_surinamus"] <- "Maschalethraupis_surinamus"
field_data$species[field_data$species == "Bucco_tamatia"] <- "Tamatia_tamatia"
field_data$species[field_data$species == "Tangara_palmarum"] <- "Thraupis_palmarum"

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

comm$t_period <- comm$Ano-1987

comm <- comm[,c(1:4,135,136,5:134)]

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
write.csv(comm, "comm_Bueno.csv")


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
# note that it excluded the islands that were sampled, but no species were captured. I will add it back manually

Arrepiado_2024 <- data.frame(
  Ilha_Ano = "Arrepiado_2024",
  Local = "Arrepiado",
  Ano = "2024",
  Data = "14/11/2024")
Arrepiado_2024[,5:118] <- 0
Arrepiado_2024$Ano <- as.integer(Arrepiado_2024$Ano)
Arrepiado_2024$Data <- as.Date(Arrepiado_2024$Data, format = "%d/%m/%Y")

Formiga_2024 <- data.frame(
  Ilha_Ano = "Formiga_2024",
  Local = "Formiga",
  Ano = "2024",
  Data = "21/06/2024")
Formiga_2024[,5:118] <- 0
Formiga_2024$Ano <- as.integer(Formiga_2024$Ano)
Formiga_2024$Data <- as.Date(Formiga_2024$Data, format = "%d/%m/%Y")

comm <- bind_rows(comm[1:7,],
  Arrepiado_2024, comm[8:26,],
  Formiga_2024, comm[27:72,])
comm <- comm[,-c(119:232)]
comm[,5:118][is.na(comm[,5:118])] <- 0

insularization_year <- "01/10/1987"
insularization_year <- as.Date(insularization_year, format = "%d/%m/%Y")

comm$t <- time_length(
  interval(insularization_year, comm$Data),
  unit = "years")

comm$t <- round(comm$t, 2)

comm$t_period <- comm$Ano-1987

comm <- comm[,c(1:4,119,120,5:118)]

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
write.csv(comm, "comm_Amarante.csv")


