########################################################
#    Mobility matrix using a radiation model
#                   PART 2
# 
# Author : Maylis Layan
# Creation date : 30/06/2020
# Last update : 30/06/2020
########################################################

# devtools::install_github('SEEG-Oxford/movement')
# library(movement)

rm(list = ls())
library(rgdal)
library(rgeos)
library(dplyr)
library(raster)
library(sf)
library(igraph)
setwd("/mnt/gaia/MMMI_Rage/Input_files/")

# Number of demes
demes <- 
#  "3demes"
  "ecoregions"

# Colors 
if (demes == "3demes") COLS = c("palegreen", "orange", "blue")
if (demes == "ecoregions") COLS = c("palegreen", "green", "darkgreen", "yellow", "orange", "lightblue", "blue")

# Aggregate the predicted mobility matrix---------------
# Load data
load("1.predict_mobility.RData")
inhab_poly = as_Spatial(st_read(paste0("gadm36_MAR_", demes, ".shp")), IDs = "AGG_ID")
if (demes == "3demes") inhab_poly = inhab_poly[order(inhab_poly$REGIONS), ]
if (demes == "ecoregions") inhab_poly = inhab_poly[order(inhab_poly$AGG_ID), ]

inhab_raster = rasterize(inhab_poly, inhab20)
plot(inhab_raster)

# Add values to outer pixels
cells = which(getValues(is.na(inhab_raster)) & getValues(!is.na(inhab20)) & getValues(inhab20 > 1))
to_consider = which(getValues(!is.na(inhab_raster)))
while (length(cells) != 0) {
  for (cell in cells) {
    adj = adjacent(inhab_raster, cell, target = to_consider)
    if (class(adj)[1] == "numeric") {
      inhab_raster[cell] = inhab_raster[adj[2]]
    } else {
      if (length(adj) != 0) {
        inhab_raster[cell] = round(mean(inhab_raster[adj[,2]], na.rm = TRUE))
      }
    }
  }
  cells = which(getValues(is.na(inhab_raster)) & getValues(!is.na(inhab20)) & getValues(inhab20 > 1))
  to_consider = which(getValues(!is.na(inhab_raster)))
}
new_poly = rasterToPolygons(inhab_raster, dissolve = TRUE)

# Correspondance between cell id and ecoregions
grouping = raster::extract(inhab20, new_poly, cellnumbers = TRUE, df = TRUE, 
                           na.rm = TRUE)

plot(inhab20)
plot(new_poly, add=T)
for (x in setdiff(grouping$cell, as.numeric(row.names(predicted_movements$movement_matrix)))) {
  rc <- rowColFromCell(inhab20, x)
  plot(extent(inhab20, rc[1], rc[1], rc[2],  rc[2]), add=TRUE, col='red', lwd=1)
}

# Mobility matrix
mobility = predicted_movements$movement_matrix
mobility = mobility[row.names(mobility) %in% grouping$cell, colnames(mobility) %in% grouping$cell]
xsum = rowsum(mobility, grouping$ID)
out = t(rowsum(t(xsum), grouping$ID))
diag(out) = 0

# Plot directed graph
# Origins corresponds to rows and destinations to columns in the movement matrix 
# See function "as.data.frame.movement_matrix" in the source code available here https://github.com/SEEG-Oxford/movement/blob/master/R/movement.R 
adjMat = matrix(1, nrow = nrow(out), ncol = ncol(out))
diag(adjMat) = 0
outLong = mutate(data.frame(out), source = row.names(out)) %>%
  tidyr::pivot_longer(-source, names_to = "destination", values_to = "weights") %>%
  mutate(destination = gsub("X", "", destination)) %>%
  filter(source != destination)
outLong$col = rep(COLS,
                  each = nrow(out)-1)

network = graph_from_adjacency_matrix(adjMat)
png(paste0("connectivity_graph_", demes, ".png"))
plot(network, vertex.color = COLS, 
     edge.color = outLong$col,
     edge.arrow.size = 0, edge.arrow.width = 0,
     edge.width = outLong$weights/90000, 
     edge.curved = 0.3, 
     layout=layout.circle)
dev.off()

# Plot polygons with same color scheme
png(paste0('connectivity_color_poly_', demes, ".png"))
plot(inhab_poly, col = COLS)
dev.off()

# Write the mobility matrix of ecoregions
write.table(out, paste0("mobility_matrix_", demes, ".txt"), sep = "\t", 
            row.names = FALSE, col.names = FALSE)
