########################################################
#    Mobility matrix using a radiation model
#                   PART 1
# 
# Author : Maylis Layan
# Creation date : 30/06/2020
# Last update : 30/06/2020
########################################################

library(rgdal)
library(raster)
library(dplyr)
#devtools::install_github('SEEG-Oxford/movement')
library(movement)

# Load raster data-------------------------------------
# Raster with number of inhabitants per pixel
inhab <- raster("MAR_ppp_v2b_2015_UNadj.tif") # available here: https://www.worldpop.org/geodata/summary?id=27690
proj4string(inhab) <- CRS("+init=epsg:4326")

# Build network on aggregated map----------------------
# Aggregation at the municipality level (20 km)
inhab20 <- raster::aggregate(inhab, 200, sum)

# Cut-off to get rid of rural areas
# Urban density = 1000 inhabitants / km2
net <- getNetwork(inhab20, min = 20000)
location_data <- data.frame(location = net$locations,
                            population = net$population,
                            x = net$coordinate[,1],
                            y = net$coordinate[,2])
location_data  <- as.location_dataframe(location_data)

# Predict movement--------------------------------------
# simulate movements (note the values of movement matrix must be integer)
predicted_flux  <- predict(originalRadiation(), location_data, symmetric = FALSE)
movement_matrix <- round(predicted_flux$movement_matrix)
movement_model <- movement(movement_matrix ~ location_data, originalRadiation())

# Predict movements
predicted_movements  <- predict(movement_model, inhab20)
save.image("../inputfiles/1.predict_mobility.RData")
