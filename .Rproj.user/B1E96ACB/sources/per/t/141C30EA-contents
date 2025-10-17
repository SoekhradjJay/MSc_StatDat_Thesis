################################################################################

# ------------------ 3. Kriging + cluster preprocessing   ----------------------

# This script is used for preprocessing the global (kriging) and local (cluster) 
# autocorrelation conditions. Make sure the relevant data and functions are 
# loaded (previous scripts).
# Last revision: 15 October 2025 (unless stated otherwise in GitHub)

################################################################################

# ---------------------------------- LIBARIES ----------------------------------

# Just make sure these packages are loaded

# Load libraries
library(sf)
library(terra)
library(sp)
library(gstat)
library(parallel)
library(data.table)

library(tidyverse)
library(dplyr)
library(raster)
library(arrow)
library(moments)

############################ Preprocessing script ##############################

# ----------- 1. create Kriging grid (called GRF grid) with id -----------------

# Make kriging grid as sf object
GRF_grid_size = 1000
GRF_grid <- st_make_grid(
  x = study_zone,
  cellsize = GRF_grid_size,
  what = "polygons",
  square = T
) %>% st_as_sf() %>%
  mutate(cell_id_GRF = 1:nrow(.))

# Also as terra object (quicker for gstat)
# as SpatialPixelsDataFrame (needed for GRF generation, later converted to sf)

# Convert to SpatialPolygons
base_grid_polygons <- as(GRF_grid, "Spatial")

# Get centroids
centroids <- coordinates(base_grid_polygons)

# Extract ID column to use in SpatialPointsDataFrame
ID_vals <- base_grid_polygons@data$cell_id_GRF

# Create SpatialPointsDataFrame with IDs
GRF_grid_sp <- SpatialPointsDataFrame(coords = centroids,
                                      data = data.frame(cell_id_GRF = ID_vals))
gridded(GRF_grid_sp) <- TRUE

# add crs to the sp object
tmp_crs = st_crs(study_zone)
proj4string(GRF_grid_sp) <- CRS(tmp_crs$proj4string)


# -------------------- 2. intersect GRF grid with target

GRF_0 <- make_GRF0(GRF_grid = GRF_grid,
                   target_zone = target_zone,
                   buffer_m = GRF_grid_size)

# ----------------- Add intersection LU and area

GRF_base <- st_intersection(GRF_0, LU_polygon)
GRF_base$area_int = st_area(GRF_base) %>% as.numeric()

# ---------------- Create a zone where clusters can be generated (near borders)

cluster_zone = generate_cluster_zone(
  original_GRF_grid = GRF_grid,
  target_zone_t = target_zone,
  bufferzone = 100
)

# -------------- Create cluster point sampling zone 
cluster_point_sampling_zone <- st_union(GRF_base)
