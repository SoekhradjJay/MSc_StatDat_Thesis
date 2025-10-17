################################################################################

# ------- 2. Load in study zone + LU map + aggregation grids  ------------------

# This script is used to load the study zone (Zuid-Holland) and create/read in 
# land use polygons. Also, aggregation grids are created here.
# Last revision: 15 October 2025 (unless stated otherwise in GitHub)

################################################################################

# ---------------------------- Load in study zone ------------------------------

# load in wijkenkaart (existing map)
prov_tmp2 = read_sf("data/study_zone/BestuurlijkeGebieden_2024.gpkg")

# make province borders (ROI/study zone)
province_border = prov_tmp2 %>%
  filter(ligt_in_provincie_naam == "Zuid-Holland") %>%
  st_make_valid() %>%
  st_union()

# make sf object
province_border = st_as_sf(province_border)

# make sure crs is 28992
st_transform(province_border, crs = 28992)

# create study_zone object 
study_zone = province_border

# --------------------------- Retrieve target zone -----------------------------

# make municipality borders
municipality_borders = prov_tmp2 %>%
  mutate(municipality_id = 1:nrow(.)) %>%
  filter(ligt_in_provincie_naam == "Zuid-Holland") %>%
  st_make_valid() %>%
  dplyr::select(municipality_id)

# make sure crs is 28992
st_transform(municipality_borders, crs = 28992)

# create target_zone object
target_zone = municipality_borders

# --------------------- Create LU_polygon (and write out) ----------------------
# choose if you want an already existing LU_polygon map or make a new one:

# number of land use polygons
n_LU_polys = 2500

# to generate a new polygon

# LU_polygon = generate_LU_polygons(study_zone, n_cells = n_LU_polys) 
# st_write(LU_polygon, "./data/LU_poly/LU_polygon.gpkg") # immediately write out

# to read a given polygon in /data/LU_poly/
LU_polygon <- st_read("data/LU_poly/LU_polygon-05082025.gpkg")

# ------------------------- Create sample grids --------------------------------

#generate grids of multiple sizes with systematic shifts (for MAUP)

plan(sequential)
step_size = 4
for(i in c(1000, 500, 100)){
  for(j in seq(0,i, i/step_size)){
    for (k in seq(0,i,i/step_size)){
      generate_grid()
      print(paste0("done with cellsize ", i, "dxdy", j, "-", k))
    }
  }
}