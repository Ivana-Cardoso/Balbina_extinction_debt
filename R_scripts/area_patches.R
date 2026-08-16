rm(list = ls())
gc()
options(scipen = 999)
set.seed(13)

# MAPBiomas Collection 9
# Load necessary packages
library(terra)
library(landscapemetrics)
library(sf)

#### 2010 ####

# Calculating the area for Marco's island using only pixel 3 MAPBiomas
# Import raster data from MAPBIOMAS 2022 Sentinel (10m)
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Rasters")
raster_2010 <- rast("mapbiomas-balbina-collection-90-2010.tif")

# Check the CRS of the imported raster
crs(raster_2010)

# Reproject the raster to UTM Zone 21S with nearest neighbor resampling
raster_2010 <- project(raster_2010, "EPSG:32721", method = "near")

# Round the raster values to avoid fractional classes
raster_2010 <- round(raster_2010)

# Check unique values in the reprojected raster
unique_values <- unique(values(raster_2010))
print(unique_values)

# Reclassify pixels 3 and 6, cause they were both forest on older versions
forest_formation_2010 <- raster_2010
forest_formation_2010[forest_formation_2010 %in% c(3, 6)] <- 3
forest_formation_2010[forest_formation_2010 != 3] <- 0

# Check if the raster has the correct properties
check_landscape(forest_formation_2010)

# If no longer needed, remove the 'Balbina' object from memory to free up space (optional)
remove(raster_2010)
plot(forest_formation_2010)

# Import environmental variables for the sites
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
sites_2010 <- read.csv("Coordenadas_UTM_Marco.csv")

# Prepare coordinates data frame from the imported sites
coordinates <- data.frame(
  X = sites_2010$X, 
  Y = sites_2010$Y, 
  id = sites_2010$ID
)

# Transform points into spatial coordinates
coordinates_wgs84 <- st_as_sf(coordinates, coords = c("X", "Y"), crs = 4326)
coordinates_utm <- st_transform(coordinates_wgs84, crs = 32721)
points(coordinates_utm, col="green", pch=21, cex = 2)

# Calculate the area of the island where each coordinate is located
area <- extract_lsm(forest_formation_2010, y = coordinates_utm, 
                    what = "lsm_p_area", extract_id = coordinates_utm$id)
area <- area[area$class == 3,]

# Rename columns to merge area and sites data frames
colnames(sites_2010)[1] = "extract_id"
sites <- merge(sites_2010, area, by = "extract_id")
sites <- sites[,c(1:5,11)]
colnames(sites)[1] <- "sites"
colnames(sites)[6] <- "area"
sites$area <- round(sites$area,2)

setdiff(sites_2010$extract_id, area$extract_id)

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
#write.csv(sites, "environment_data_Marco2010.csv")


#### 2015 ####
rm(list = ls())
gc()

# Calculating the area for Marco's island using only pixel 3 MAPBiomas
# Import raster data from MAPBIOMAS 2022 Sentinel (10m)
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Rasters")
raster_2015 <- rast("mapbiomas-balbina-collection-90-2015.tif")

# Check the CRS of the imported raster
crs(raster_2015)

# Reproject the raster to UTM Zone 21S with nearest neighbor resampling
raster_2015 <- project(raster_2015, "EPSG:32721", method = "near")

# Round the raster values to avoid fractional classes
raster_2015 <- round(raster_2015)

# Check unique values in the reprojected raster
unique_values <- unique(values(raster_2015))
print(unique_values)

# Reclassify pixels 3 and 6, cause they were both forest on older versions
forest_formation_2015 <- raster_2015
forest_formation_2015[forest_formation_2015 %in% c(3, 6)] <- 3
forest_formation_2015[forest_formation_2015 != 3] <- 0

# Check if the raster has the correct properties
check_landscape(forest_formation_2015)

# If no longer needed, remove the 'Balbina' object from memory to free up space (optional)
remove(raster_2015)
plot(forest_formation_2015)

# Import environmental variables for the sites
setwd("C:/Users/ivana/OneDrive/PhD_INPA/Dados_Balbina/KNB_Anderson")
sites_2015 <- read.csv("balbina_environmental_variables.csv")

# Prepare coordinates data frame from the imported sites
coordinates <- data.frame(
  X = sites_2015$longitude.WGS84, 
  Y = sites_2015$latitude.WGS84, 
  id = sites_2015$site
)

# Transform points into spatial coordinates
coordinates_wgs84 <- st_as_sf(coordinates, coords = c("X", "Y"), crs = 4326)
coordinates_utm <- st_transform(coordinates_wgs84, crs = 32721)
points(coordinates_utm, col = "red")

# Calculate the area of the island where each coordinate is located
area <- extract_lsm(forest_formation_2015, y = coordinates_utm, 
                    what = "lsm_p_area", extract_id = coordinates_utm$id)
area <- area[area$class == 3,]

# Rename columns to merge area and sites data frames
colnames(sites_2015)[1] = "extract_id"
sites <- merge(sites_2015, area, by = "extract_id")
sites <- sites[,c(1:4,50)]
colnames(sites)[1] <- "sites"
colnames(sites)[5] <- "area"

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
#write.csv(sites, "environment_data_Anderson2015.csv")


#### 2023 ####
rm(list = ls())
gc()

# Calculating the area for Marco's island using only pixel 3 MAPBiomas
# Import raster data from MAPBIOMAS 2022 Sentinel (10m)
setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Rasters")
raster_2023 <- rast("00mapbiomas-balbina-collection-90-2023.tif")

# Check the CRS of the imported raster
crs(raster_2023)

# Reproject the raster to UTM Zone 21S with nearest neighbor resampling
raster_2023 <- project(raster_2023, "EPSG:32721", method = "near")

# Round the raster values to avoid fractional classes
raster_2023 <- round(raster_2023)

# Check unique values in the reprojected raster
unique_values <- unique(values(raster_2023))
print(unique_values)

# Reclassify pixels 3 and 6, cause they were both forest on older versions
forest_formation_2023 <- raster_2023
forest_formation_2023[forest_formation_2023 %in% c(3, 6)] <- 3
forest_formation_2023[forest_formation_2023 != 3] <- 0

# Check if the raster has the correct properties
check_landscape(forest_formation_2023)

# If no longer needed, remove the 'Balbina' object from memory to free up space (optional)
remove(raster_2023)
plot(forest_formation_2023)

# Import environmental variables for the sites
setwd("C:/Users/ivana/OneDrive/PhD_INPA/1.Diversity_question/Data")
sites_2023 <- read.csv("Amarante_environmental_variables.csv")

# Prepare coordinates data frame from the imported sites
coordinates <- data.frame(
  X = sites_2023$longitude.WGS84, 
  Y = sites_2023$latitude.WGS84, 
  id = sites_2023$site
)

# Transform points into spatial coordinates
coordinates_wgs84 <- st_as_sf(coordinates, coords = c("X", "Y"), crs = 4326)
coordinates_utm <- st_transform(coordinates_wgs84, crs = 32721)
points(coordinates_utm, col = "red")

# Calculate the area of the island where each coordinate is located
area <- extract_lsm(forest_formation_2023, y = coordinates_utm, 
                    what = "lsm_p_area", extract_id = coordinates_utm$id)
area <- area[area$class == 3,]

# Rename columns to merge area and sites data frames
colnames(sites_2023)[1] = "extract_id"
sites <- merge(sites_2023, area, by = "extract_id")
sites <- sites[,c(1:4,50)]
colnames(sites)[1] <- "sites"
colnames(sites)[5] <- "area"

setwd("C:/Users/ivana/OneDrive/PhD_INPA/3.Extinction_debt/Dados_Marco")
#write.csv(sites, "environment_data_Amarante2023.csv")
