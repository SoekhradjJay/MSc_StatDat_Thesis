################################################################################

# ------- 4. Run simulations   ------------------

# This script is used to run the simulations.
# Last revision: 15 October 2025 (unless stated otherwise in GitHub)

################################################################################

library(sf)
library(terra)
library(sp)
library(gstat)
library(future)
library(future.apply)
library(parallel)
library(data.table)

library(tidyverse)
library(dplyr)
library(raster)
library(arrow)
library(moments)

################################################################################

# -------------------- Simulation conditions  ----------------------------------

# create 2x2x2 simulation conditions. Each row is a scenario.
sim_conditions <- expand.grid(
  c("smooth", "nonsmooth"), # Kriging (on/off)
  c("informative", "noninformative"), # land use weights (strong vs. flat)
  c("clustered", "nonclustered") # cluster presence
)

# rename the column names
names(sim_conditions) <- c("GRF", "LU_info", "clustering")

# get cores
workers = detectCores() - 1

################################################################################

# ------ first step create in parallel n_block GRFs (kriging fields)

# make sure we are in non-parallel mode
plan(sequential)

# start timer
start = Sys.time()

# reps per condition (number of point datasets per condition)
n_block = 100

# when no kriging is performed, all weights are 1:
# already create GRF_list when no kriging is performed
grf_list_noGRF <- lapply(1:n_block, function(i)
  rep(1, GRF_grid_sp@data %>% nrow()))

# get frame ready to store the rRMSE and other metrics
all_results = tibble()


# perform simulation for all conditions
for (condition_i in 1:8) {
  
  print(paste0("start condition,", condition_i))
  
  # gather simulation conditions for each scenario
  current_condition = condition_i
  is_smooth = sim_conditions[current_condition, 1]
  is_informative = sim_conditions[current_condition, 2]
  is_clustered = sim_conditions[current_condition, 3]
  
  # ------------------------- CONDITION 1: KRIGING ------------------------- 
  
  # set to non-parallel within the loop
  plan(sequential) 
  
  # Perform kriging (or not based on is_smooth)
  if (is_smooth == "smooth"){
    plan(multisession, workers = workers)
    grf_list <- future_lapply(
      1:n_block,
      FUN = function(i) {
        simulate_grf(
          base_grid_sp = GRF_grid_sp,
          study_zone_s = study_zone,
          is_smooth = is_smooth
        )
      },
      future.seed = 7
    )
  } else {
    grf_list = grf_list_noGRF
  }

  # set to non-parallel within the loop
  plan(sequential)
  
  # get all kriging weights into a matrix (geometries are constant)
  GRF_mat <- make_GRF_mat(
    GRFs = grf_list,
    n_block = n_block,
    GRF_grid_sp_ncell = GRF_grid_sp@data$cell_id_GRF %>% length()
  )
  
  # -------------------------- running point sim + SAW/DASYM
  
  # tibble to store only one run of rRMSE
  output_RMSE = tibble()
  
  # Get the land use weights
  # in next version: can be put in preprocessing
  x_vals = seq(0.01, 0.5, length.out = 8) # how many classes are there?
  alphas = c(0, 0.1, 0.2, 0.5, 0.7, 1, 2, 5, 7, 10) # how many lambda sets?
  all_trial_lambdas = lapply(alphas, generate_lambda_LU_set, x_vals = x_vals)
  
  # ------------------- CONDITION 2/3: LAND USE + CLUSTERS -------------------
  # note that for all repetitions within a scenario, the kriging grids have 
  # already been calculated. Now, in parallel, we simulate the land use and 
  # cluster effects.
  
  for (i in 1:n_block) {
    print(paste0("start n_block, ", i))
    
    # first add land use to the dataset
    GRF_df <- preprocess_GRF(
      i = i, 
      GRF_mat = GRF_mat, 
      GRF_base = GRF_base, # from preprocessing
      informative = is_informative,
      nrow_LU_polys = n_LU_polys,
      n_LU_classes = 8,
      lambda_sets_obj = all_trial_lambdas
    )
    
    # calculate density per intersection and 
    # convert to normalized baseline count
    GRF_df <- calculate_npop(
      GRF_df_full = GRF_df,
      n_cluster = 400,
      size_cluster = 50,
      clustered = is_clustered,
      pop_total = 200000
    )
    
    # ------------------- point generation:
    # creating a point dataset takes 1.5 minutes appr.
    print("start making points")
    
    # set-up parallel environment
    plan(multisession, workers = workers)
    
    # holder for all point geometries
    samples_list = list()
    
    # Sample points in parallel and tag with original 
    # row index (row = polygon geometry)
    
    samples_list <- future_lapply(1:nrow(GRF_df), function(j) {
      n <- GRF_df$norm_npop_unbiased[j]
      if (n > 0) {
        st_sample(GRF_df[j, ], size = n)
      } else {
        NULL
      }
    }, future.seed = 7)
    
    # set back to non-parallel mode
    plan(sequential)
    
    # add all point geometries to list
    samples_list = samples_list[!sapply(samples_list, is.null)]
    samples = do.call(c, samples_list) %>% st_sf() 
    
    # add a cluster id (not necessary, but for overview)
    samples$belongs_to_cluster = 0 # these not
    samples$cluster_id = NA # those with 0, do not have an ID
    
    # generate clusters
    if (is_clustered == "clustered") {
      all_cluster_points = generate_clusters(
        cluster_zone = cluster_zone,
        CL_sampling_zone = cluster_point_sampling_zone,
        n_cluster = 400,
        size_cluster = 50
      ) %>% rename(geometry = x)
      all_sampled_points = rbind(all_cluster_points, samples)
    } else{
      all_sampled_points = samples
    }
    
    # point generation completed
    print("done with points, start with SAW/DASYM")
    
    # ------------------- point aggregation + SAW/DASYM
    # In the preprocessing, we generated aggregation grids. We open these 
    # grids for point aggregation
    
    
    # Open the aggregation grids
    ds <- open_dataset("data/grids",
                       partitioning = c("cell_size", "dx_id", "dy_id"))
    
    # Extract available tiles (one row per file/partition)
    tiles <- ds %>%
      dplyr::select(cell_size, dx_id, dy_id) %>%
      distinct() %>%
      collect() %>% filter(dx_id != cell_size &
                             dy_id != cell_size) %>%
      arrange(cell_size, dx_id, dy_id)
    
    # create a frame were LU_id and LU_class are together (from GRF_df)
    LU_poly_filled =
      GRF_df %>% dplyr::select(LU_id, LU_class)
    
    # prepare parallel processing
    plan(multisession, workers = workers)
    
    # Split the job into multiple tiles for parallel processing
    # Each tile is a grid
    tile_list <- split(tiles, seq_len(nrow(tiles)))  
    
    # Perform point aggregation and immediately perform SAW and DASYM
    results <- future_lapply(tile_list, function(tile_row) {
      process_tile(tile_row = tile_row,
                   condition_id = condition_i,
                   rep_id = i, 
                   LU_polys = LU_poly_filled, 
                   points = all_sampled_points)
    }, future.seed = 7)
    
    # Reset to non-parallel when done
    plan(sequential)  
    
    output_RMSE = bind_rows(output_RMSE, bind_rows(results))
  }
  
  # write output of one condition (x many blocks) into a file
  file_path = file.path("data/results/", sprintf("condition_%03d.parquet", condition_i))
  write_parquet(output_RMSE, file_path)
  
  # append results to an all_results frame (in user session)
  #all_results = bind_rows(all_results, output_RMSE)
}
end = Sys.time()
end - start

