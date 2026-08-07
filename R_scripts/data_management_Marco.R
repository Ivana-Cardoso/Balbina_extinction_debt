rm(list = ls())
gc()
options(scipen = 999)
set.seed(13)

# Set working directory
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco/Dados_originais")

# Import bird capture data
field_data <- read.csv("capturas1.csv", header = TRUE)

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

# Create occurrence column
field_data$occurrence <- 1

# Construct a site-by-species matrix for island sites
comm <- as.data.frame(
  tapply(
    field_data$occurrence,
    list(field_data$Local, field_data$Espécie),
    sum
  )
)

# Replace NA values with zero
comm[is.na(comm)] <- 0

# Sum the transects that are on the same island
comm[7, ]  <- colSums(comm[c(7,8), ])
comm[16, ] <- colSums(comm[c(16,17), ])
comm[19, ] <- colSums(comm[c(19,20), ])
comm[24, ] <- colSums(comm[c(24,25), ])

# Remove the second line of each pair
comm <- comm[-c(8, 17, 20, 25),]

# Rename the lines that were summed
row.names(comm)[7] <- "1815"
row.names(comm)[9] <- "213"
row.names(comm)[15] <- "475"
row.names(comm)[17] <- "690"
row.names(comm)[20] <- "FCA_LOR"
row.names(comm)[21] <- "FCC_WB"

# Remove species name without epiteto especifico
comm <- comm[, -c(46,55,62, 89)]

# Some of Marco's islands were sampled by Anderson and me, so I will add a new 
# column with our islands names
comm$island = NA
comm = comm[,c(96,1:95)]
comm$island[9] = "Cipoal"
comm$island[3] = "Piquia"
comm$island[6] = "Coata"
comm$island[15] = "Martelo"
comm$island[7] = "Gaviao_real"
comm$island[4] = "Abusado"
comm$island[14] = "Panema"
comm$island[17] = "Beco_do_catitu"
comm$island[2] = "Pontal"
comm$island[18] = "Relogio"
comm$island[19] = "Sapupara"
comm$island[20] = "CF_Loreno"
comm$island[21] = "CF_WABA"

# Marco used the area of the island as IDs
comm$area = rownames(comm)
comm = comm[, c(1,97,2:96)]

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
write.csv(comm, "comm_Marco.csv")
