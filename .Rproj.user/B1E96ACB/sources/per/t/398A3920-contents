################################################################################

# --------------------------- 1. Load functions  -------------------------------

# This script is used to load all functions/formulas for the simulation script
# Last revision: 15 October 2025 (unless stated otherwise in GitHub)

################################################################################

# ---------------------------------- LIBARIES ----------------------------------

# Just make sure these packages are loaded

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

# ---------------------------------- FORMULAS --------------------------------


# ***************************** sim preparation *******************************

# to generate the land use polygon geometries with a LU_id (identifier)
generate_LU_polygons <- function(big_polygon, n_cells) {
  # 1. Sample random points inside the polygon
  sample_points <- st_sample(big_polygon, size = n_cells, type = "random")
  
  # 2. Generate Voronoi polygons around the points
  envelope <- st_as_sfc(st_bbox(big_polygon))
  voronoi_raw <- st_voronoi(st_union(sample_points), envelope = envelope)
  voronoi_polys <- st_collection_extract(voronoi_raw, "POLYGON")
  voronoi_sf <- st_sf(geometry = voronoi_polys)
  
  # 3. Clip Voronoi cells to the main polygon
  voronoi_clipped <- st_intersection(voronoi_sf, big_polygon)
  
  # under revision (could be added)
  # 4. Optionally apply concaveman to random subset for irregular shapes (under development)
  # apply_concave <- sample(c(TRUE, FALSE), size = nrow(voronoi_clipped), replace = TRUE)
  # voronoi_clipped$geometry <- map2(voronoi_clipped$geometry, apply_concave, function(geom, use_concave) {
  #   if (use_concave && length(st_coordinates(geom)) >= 3) {
  #     tryCatch({
  #       concaveman(geom)
  #     }, error = function(e) {
  #       geom  # fallback to original
  #     })
  #   } else {
  #     geom
  #   }
  # }) %>% st_sfc(crs = st_crs(big_polygon))
  
  # 5. Add ID
  voronoi_clipped$LU_id = 1:nrow(voronoi_clipped)
  
  # Return diverse subpolygons
  st_sf(voronoi_clipped)
}

# to generate grids over which we aggregate points. It is meant as function to
# pass to parallel processing, meaning that cell_size, dx (shift in x), dy
# (shift in y) can be looped over.

generate_grid = function(cellsize = i,
                         dx = j,
                         dy = k,
                         study_zone_s = study_zone,
                         target_zone_t = target_zone) {
  # 1. Define grid parameters
  crs_rd <- "EPSG:28992"
  
  # 2. Convert bounding area to SpatVector (terra-compatible)
  province_vect <- vect(study_zone_s)  # convert sf -> terra
  ext <- ext(province_vect)  # get extent
  
  # create bigger bounding box (bottom row and left column)
  ext[1] <- as.numeric(ext[1]) - cellsize
  ext[3] <- as.numeric(ext[3]) - cellsize
  ext[2] <- as.numeric(ext[2]) + cellsize
  ext[4] <- as.numeric(ext[4]) + cellsize
  
  # 3. Create raster template
  r <- rast(
    extent = ext,
    resolution = cellsize,
    crs = crs(province_vect)
  )
  
  # add origin-shift
  origin(r)[1] <- origin(r)[1] + dx # change origin
  origin(r)[2] <- origin(r)[2] + dy # change origin
  
  # 4. Convert raster cells to polygons (grid)
  shifted_grid <- as.polygons(r) |> st_as_sf()
  shifted_grid$cell_id <- 1:nrow(shifted_grid)
  
  # Filter only relevant cells (intersected cells)
  
  # Setup parallel plan to do intersections grid and municipality
  plan(multisession, workers = detectCores() -1)
  
  intersection_list <- future_lapply(1:nrow(target_zone_t), function(i) {
    muni_i <- target_zone_t[i, ]
    inter <- st_intersection(shifted_grid, muni_i)
    
    if (nrow(inter) > 0) {
      inter$municipality_id <- muni_i$municipality_id
    }
    
    inter
  }, future.seed = 7)
  plan(sequential)
  
  # Retrieve intersected cells
  attr_list <- lapply(intersection_list, st_set_geometry, NULL) # get attributes
  ids_only = do.call(rbind, attr_list) # get ids
  grid_freq = table(ids_only$cell_id) # count the cell_id
  single_ids <- names(grid_freq[grid_freq == 1]) # which ones have exactly 1?
  ids_only <- ids_only[!(ids_only$cell_id %in% single_ids), ] # filter those out
  
  # filter only relevant cells from the entire grid
  relevant_cells = shifted_grid[shifted_grid$cell_id %in% ids_only$cell_id, ]
  
  # intersect those cells again (faster than do.call on intersection_list).
  intersection_st = st_intersection(relevant_cells, target_zone_t)
  
  # add area to intersection source and target
  intersection_st$area_st = st_area(intersection_st) %>% as.numeric()
  
  # add an id to each intersection
  intersection_st$st_id = 1:nrow(intersection_st)
  
  ############################ add metadata
  intersection_st = intersection_st %>%
    mutate(cell_size = cellsize,
           dx_id = dx,
           dy_id = dy)
  
  ########################## prepare for write_out
  intersection_st$wkt <- st_as_text(intersection_st$geometry)
  intersection_st_out <- intersection_st %>% dplyr::select(-geometry) %>% st_drop_geometry()
  
  # write out
  write_dataset(
    intersection_st_out,
    path = "./data/grids/",
    format = "parquet",
    partitioning = c("cell_size", "dx_id", "dy_id"),
    compression = "zstd"
  )
}

# ***************************** sim preprocessing *****************************

# to intersect the Kriging grid called GRF_grid (invoked in separate script) 
# with the target_zone and create a buffer of buffer_m around 
# the cells that are intersected
make_GRF0 <- function(GRF_grid = GRF_grid,
                      target_zone = target_zone,
                      buffer_m) {
  ###### First find cells intersected by target borders
  # intersect GRF grid with border
  GRF_0_prep = st_intersection(GRF_grid, target_zone)
  
  # get back which cells are actually intersected
  grid_freq = table(GRF_0_prep$cell_id_GRF)
  single_ids <- names(grid_freq[grid_freq == 1])
  is_intersected_munborder <- GRF_0_prep[!(GRF_0_prep$cell_id_GRF %in% single_ids), ]
  
  ##### Then, create buffer zone around intersected cells
  
  # get buffer distance (choose wisely (1x GRF grid cell size recommended))
  buffer_distance_m = buffer_m
  tmp_buf =  st_buffer(is_intersected_munborder, dist = buffer_distance_m)
  
  # join buffer into one object
  merged_buffer <- st_union(tmp_buf)
  
  # Find intersections buffer and original grid
  ints <- st_intersects(merged_buffer, GRF_grid, sparse = TRUE)
  
  # Find the indices (rows) of grid cells that intersect any buffered cell
  intersected_idx <- unique(unlist(ints))
  
  # Subset the already intersected grid on future intersected cells + buffered cells
  GRF_0_whole_cells <- GRF_grid[intersected_idx, ]
  
  # Now re-intersect the filtered grid with the target zone
  GRF_0 <- st_intersection(GRF_0_whole_cells, target_zone)
  
  return(GRF_0)
}

# to create a zone where clusters may appear, which is around intersected cells
generate_cluster_zone <- function(original_GRF_grid = GRF_grid,
                                  target_zone_t = target_zone,
                                  bufferzone = 100) {
  # 1. Select epicenters of cluster locations within buffer zone
  
  # intersect GRF grid with border
  GRF_0_prep = st_intersection(original_GRF_grid, target_zone_t)
  
  # get back which cells are actually intersected
  grid_freq = table(GRF_0_prep$cell_id_GRF)
  single_ids <- names(grid_freq[grid_freq == 1])
  is_intersected_munborder <- GRF_0_prep[!(GRF_0_prep$cell_id_GRF %in% single_ids), ]
  
  # get the union around is_intersected to sample points close to border
  tmp_buf_simpoints =  st_buffer(is_intersected_munborder, dist = bufferzone) %>% st_union()
}

# **************************** running simulation ****************************

# to get the kriging weight called lambda_GRF as vector out. 
# Each index corresponds to the cell indicated by the GRF_grid_sp cell order. 
# This is meant to be run the simulation loop.

simulate_grf <- function(base_grid_sp = GRF_grid_sp,
                         study_zone_s = study_zone,
                         is_smooth) {
  # Step 1: Sample random points in study zone
  pts <- st_sample(study_zone_s, size = 1000)
  pts <- st_as_sf(pts)
  pts$z <- rnorm(nrow(pts))
  
  # Step 2: Convert to SpatialPointsDataFrame
  pts_sp <- as(pts, "Spatial")
  
  # Step 3: Variogram (exponential model)
  range_val <- if (is_smooth == "smooth")
    5000
  else
    500
  vgm_model <- vgm(
    psill = 1, # Q
    model = "Exp",
    range = range_val,
    nugget = 0.1 # Q_0
  )
  
  # Step 4: Execute Kriging
  g_model <- gstat(formula = z ~ 1,
                   data = pts_sp,
                   model = vgm_model)
  kriged <- predict(g_model, newdata = base_grid_sp)
  
  # Step 5: Return predicted values
  out = kriged$var1.pred
  return(out) # these are the weights
}

# to put the output of simulate_grf in a matrix style with indicator which grid
# cell and which simulation block lambda_GRF belongs to.
make_GRF_mat <- function(GRFs, n_block, GRF_grid_sp_ncell) {
  ### get the columns ready
  
  # create one big column with GRFs
  lambda_GRF <- do.call(c, GRFs)
  cell_id_GRF <- rep(1:GRF_grid_sp_ncell, n_block)
  GRF_ind <- rep(1:n_block, each = GRF_grid_sp_ncell)
  
  # prepare matrix
  GRF_matrix = matrix(nrow = length(lambda_GRF), ncol = 3)
  
  # fill matrix
  GRF_matrix[, 1] <- GRF_ind
  GRF_matrix[, 2] <- cell_id_GRF
  GRF_matrix[, 3] <- lambda_GRF
  
  # give column names
  colnames(GRF_matrix) <- c("GRF_ind", "cell_id_GRF", "lambda_GRF")
  
  #output
  
  return(GRF_matrix)
}

# to give each LU_polygon a random LU_class and weight. ouputs this in a df
# ready for left_joining to LU_polygon. Used in wrapper preprocess_GRF.
create_lambdaLU <- function(LUinfo, nrow_LU_polys = n_LU_polys, n_LU_classes, lambda_sets_obj = all_trial_lambdas) {
  # create output df with already random LU_classes per LU polygon
  LU_df = data.frame(
    LU_id = 1:nrow_LU_polys,
    LU_class = sample(1:n_LU_classes, nrow_LU_polys, replace = T)
  )
  
  # join necessary lambda_LU based on informativeness LU
  if (LUinfo == "informative") {
    lambda_LU_informative_df = data.frame(
      LU_class = 1:n_LU_classes,
      lambda_LU = lambda_sets_obj[[length(lambda_sets_obj)]]
    )
    out = left_join(LU_df, lambda_LU_informative_df, by = "LU_class")
  } else {
    lambda_LU_noninformative = data.frame(
      LU_class = 1:n_LU_classes,
      lambda_LU = c(1, 1, 1, 1, 1, 1, 1, 1)
    )
    out = left_join(LU_df, lambda_LU_noninformative, by = "LU_class")
  }
  return(out)
}

# wrapper function. Takes as input make_GRF_mat output and adds to GRF_base
# data frame made
preprocess_GRF <- function(i,
                           GRF_mat,
                           GRF_base,
                           informative,
                           nrow_LU_polys,
                           n_LU_classes,
                           lambda_sets_obj = all_trial_lambdas) {
  # selection of the right GRF
  GRF_selected <- GRF_mat[GRF_mat[, 1] == i, ]
  
  # add GRF to the base frame (geometries)
  GRF_df <- left_join(GRF_base, as.data.frame(GRF_selected), by = "cell_id_GRF")
  
  # create a LU_class configuration
  lambdaLU_df <- create_lambdaLU(LUinfo = informative,
                                 # depends on condition
                                 nrow_LU_polys = nrow_LU_polys,
                                 n_LU_classes = n_LU_classes, 
                                 lambda_sets_obj = all_trial_lambdas)
  
  # add to the GRF_df
  GRF_df <- left_join(GRF_df, lambdaLU_df, by = "LU_id")
  
  # return preprocessed frame
  
  return(GRF_df)
}

# to calculate number of people per intersection GRF grid and land use poly
# takes output from preprocess_GRF as GRF_df_full argument.
calculate_npop <- function(GRF_df_full,
                           n_cluster = 20,
                           size_cluster = 100,
                           clustered,
                           pop_total = 1000000) {
  # 1. Calculate the density per polygon
  GRF_df = GRF_df_full %>%
    mutate(base_density = exp(lambda_GRF) * lambda_LU) %>%
    mutate(npop = base_density * area_int)
  
  # 2. normalize the population to
  totpop = sum(GRF_df$npop) # the actual number from simulation
  
  if (clustered == "clustered") {
    pop_total = pop_total - (n_cluster * size_cluster)
  } else{
    pop_total = pop_total
  }
  
  GRF_df = GRF_df %>%
    mutate(norm_npop = npop / totpop * pop_total)
  
  # unbiased rounding (for totals)
  x <- GRF_df$norm_npop
  # Floor all values
  x_floor <- floor(x)
  
  # Compute how many units are left to distribute
  remainder <- sum(x) - sum(x_floor)
  to_add <- round(remainder)
  
  # Compute the decimal parts
  decimals <- x - x_floor
  
  # Order indices by largest decimal parts
  order_idx <- order(decimals, decreasing = TRUE)
  
  # Add 1 to the top 'to_add' indices
  x_rounded <- x_floor
  if (to_add > 0) {
    x_rounded[order_idx[1:to_add]] <- x_rounded[order_idx[1:to_add]] + 1
  }
  
  # the unbiased rounding in frame
  GRF_df$norm_npop_unbiased = x_rounded
  
  # return GRF_df
  return(GRF_df)
}

# to calculate contribution of lambda_LU and lambda_GRF to count
generate_lambda_LU_set <- function(alpha, x_vals){
  lambda = exp(-alpha * x_vals)
  norm_lambda = lambda/sum(lambda, na.rm = T)
  return(norm_lambda)
}


# to generate the real point geometries. this does not add the cluster yet.
generate_base_points <- function(GRF_df, clustered) {
  # Calculate the points (as baseline)
  samples_list <- lapply(1:nrow(GRF_df), function(i) {
    print(i)
    n <- round(GRF_df$norm_npop[i], 0)
    if (n > 0)
      st_sample(GRF_df[i, ], size = n)
    else
      NULL
  })
  
  return(samples)
}

# to generate the cluster point geometries.
generate_clusters <- function(cluster_zone = cluster_zone,
                              CL_sampling_zone = cluster_point_sampling_zone,
                              n_cluster,
                              size_cluster) {
  # 1. Sample epicenters of cluster locations within buffer zone
  
  # sample some epicenters
  epis <- st_sample(cluster_point_sampling_zone, n_cluster)
  
  # 2. Perform actual clustering with epicenters
  
  # for rejection sampling, get domain for points to fall in
  tmp_union = CL_sampling_zone
  
  all_cluster_points <- lapply(1:length(epis), function(i) {
    oversampling_factor <- 2  # initial value
    max_tries <- 10           # safety stop to avoid infinite loop
    try_count <- 0
    
    repeat {
      try_count <- try_count + 1
      
      # 1. Simulate oversampled points
      n_points <- size_cluster * oversampling_factor
      clust_coordinates <- mvtnorm::rmvnorm(
        n = n_points,
        mean = st_coordinates(epis[i]) %>% as.numeric(),
        sigma = matrix(c(2000, 0, 0, 2000), nrow = 2)
      )
      colnames(clust_coordinates) <- c("X", "Y")
      
      # 2. Convert to sf
      clust_sf <- st_as_sf(
        as.data.frame(clust_coordinates),
        coords = c("X", "Y"),
        crs = 28992,
        sf_column_name = "x"
      )
      
      # 3. Rejection: keep only inside tmp_union
      clust_sf <- clust_sf[st_within(clust_sf, tmp_union, sparse = FALSE), ]
      
      # 4. Check if we have enough points
      if (nrow(clust_sf) >= size_cluster) {
        clust_sf <- clust_sf[sample(nrow(clust_sf), size_cluster), ]
        break
      } else {
        oversampling_factor <- ceiling(oversampling_factor * 1.5)
        if (try_count >= max_tries) {
          warning(
            paste(
              "Cluster",
              i,
              "did not reach",
              size_cluster,
              "points after",
              try_count,
              "tries. Returning",
              nrow(clust_sf),
              "points."
            )
          )
          break  # fallback: return whatever we have
        }
      }
    }
    
    clust_sf$cluster_id <- i
    return(clust_sf)
  })
  
  # 3. Prepare the output
  # get all cluster points together
  all_cluster_points <- do.call(rbind, all_cluster_points)
  
  # add columns to distinguish between cluster/noncluster points
  all_cluster_points$belongs_to_cluster = 1 # these belong to cluster
  
  return(all_cluster_points)
}

# ************************* Apply SAW/DASYM and get RMSE ***********************

# to generate lambdas to test in the simulation for DASYM
# deprecated, but could be useful as 'other way' to generate land use weights
generate_similar_lambda <- function(lambda_ref, noise_strength = 0.1) {
  lambda_ref <- lambda_ref / sum(lambda_ref)
  noise <- runif(length(lambda_ref))
  noise <- noise / sum(noise)
  
  lambda_sim <- (1 - noise_strength) * lambda_ref + noise_strength * noise
  return(lambda_sim / sum(lambda_sim))
}

# to perform SAW and DASYM but as a split tile, meaning in parallel process.
# meant for used in a loop.
process_tile <- function(tile_row,
                         LU_polys = LU_poly_filled,
                         points = all_sampled_points,
                         lambdas = all_trial_lambdas,
                         condition_id,
                         rep_id) {
  ########################## Get a grid based on tile_row
  ds <- open_dataset("data/grids", partitioning = c("cell_size", "dx_id", "dy_id"))
  df <- ds %>%
    collect() %>%
    filter(cell_size == tile_row$cell_size,
           dx_id == tile_row$dx_id,
           dy_id == tile_row$dy_id)
  
  ########################### Preprocessing (point aggregation)
  
  #------------ convert back to SF
  df = df %>%
    mutate(geometry = st_as_sfc(wkt)) %>%
    st_as_sf(crs = 28992) %>%
    dplyr::select(-wkt)
  
  # ----------- intersect with land use polygons + add area and id
  df_full = st_intersection(df, LU_polys)
  df_full$stc_id = 1:nrow(df_full)
  df_full$area_stc = st_area(df_full) %>% as.numeric()
  
  # ----------- aggregate counts to y_stc + create final frame
  tmp_join = st_join(x = st_as_sf(all_sampled_points),
                     y = df_full,
                     join = st_within)
  
  final_frame =  table(tmp_join$stc_id) %>%
    as.data.frame() %>%
    rename(stc_id = Var1, y_stc_true = Freq) %>%
    left_join(
      tmp_join %>%
        st_drop_geometry() %>%
        mutate(stc_id = as.factor(stc_id)) %>%
        distinct(stc_id, .keep_all = TRUE),
      by = "stc_id"
    ) %>% as.data.table()
  
  ########################## Perform SAW
  
  # calculate necessary ingredients for SAW
  area_s <- final_frame[, .(area_s = sum(area_stc)), by = cell_id]
  area_st <- final_frame[, .(area_st_1 = sum(area_stc)), by = st_id]
  y_s     <- final_frame[, .(y_s = sum(y_stc_true)), by = cell_id]
  y_t_true <- final_frame[, .(y_t_true = sum(y_stc_true)), by = municipality_id]
  
  # add everything to the frame
  tmp1 <- left_join(final_frame, area_s, by = "cell_id")
  tmp2 <- left_join(tmp1, area_st, by = "st_id")
  tmp3 <- left_join(tmp2, y_s, by = "cell_id")
  
  # estimation y_st
  tmp4 = tmp3 %>%
    mutate(y_st_estim = y_s * (area_st_1 / area_s))
  
  # gather error into frame
  error_frame_SAW = tmp4[, .(y_t_estim = sum(y_st_estim)), by = municipality_id]
  final_error_frame_SAW = left_join(error_frame_SAW, y_t_true, by = "municipality_id")
  
  # prepare output
  estims_SAW <- evaluate_estimates(
    df = final_error_frame_SAW,
    estimate_col = "y_t_estim",
    truth_col = "y_t_true",
    method = paste0("SAW"),
    dx = tile_row$dx_id,
    dy = tile_row$dy_id,
    cell_size = tile_row$cell_size
  )
  
  ########################## Perform DASYM
  
  # calculate necessary ingredients for DASYM
  y_s     <- final_frame[, .(y_s = sum(y_stc_true)), by = cell_id]
  y_t_true <- final_frame[, .(y_t_true = sum(y_stc_true)), by = municipality_id]
  area_sc <- final_frame[, .(area_sc = sum(area_stc)), by = .(cell_id, LU_class)]
  
  # repeat this for each lambda_trial
  output_dasym = list()
  for (i in 1:length(all_trial_lambdas)) {
    reference_lambdas = all_trial_lambdas[[i]]
    
    # imagine you have an external output for the lambda vec
    lambda_df <- data.frame(LU_class = 1:8, trial_lambda = reference_lambdas)
    
    area_sc_full <- left_join(area_sc, lambda_df, by = "LU_class")
    weight_frame <- area_sc_full[, .(weight_denom = sum(trial_lambda * area_sc)), by = cell_id]
    
    # create dasym frame
    final_dasym_frame_tmp <- left_join(final_frame, weight_frame, by = "cell_id") # add w_c (dasym weights)
    final_dasym_frame = left_join(final_dasym_frame_tmp, lambda_df, by = "LU_class") # get the trial_lambdas
    final_dasym_frame_final = left_join(final_dasym_frame, y_s, by = "cell_id") # add total for source

    # estimation y_stc
    tmp_dasym_out = final_dasym_frame_final %>%
      mutate(y_stc_estim = y_s * area_stc * trial_lambda / weight_denom) 
    # make error frame
    error_frame_dasym = tmp_dasym_out[, .(y_t_estim = sum(y_stc_estim)), by = municipality_id]
    final_error_frame_dasym = left_join(error_frame_dasym, y_t_true, by = "municipality_id")
    
    #apply estimates
    estims_dasym <- evaluate_estimates(
      df = final_error_frame_dasym,
      estimate_col = "y_t_estim",
      truth_col = "y_t_true",
      method = paste0("dasym_", i),
      dx = tile_row$dx_id,
      dy = tile_row$dy_id,
      cell_size = tile_row$cell_size
    )
    
    
    # put the looped vecs in list
    output_dasym[[i]] <- estims_dasym
  }
  
  # prepare_output
  out <- rbind(estims_SAW, bind_rows(output_dasym))
  out$condition_id = condition_id
  out$rep_id = rep_id
  
  return(out)
}

# to calculate the RMSE and more metrics
evaluate_estimates <- function(df,
                               estimate_col,
                               truth_col,
                               method,
                               dx,
                               dy,
                               cell_size) {
  # Extract vectors
  y_hat <- df[[estimate_col]]
  y_true <- df[[truth_col]]
  
  # Residuals
  error <- y_hat - y_true
  
  tibble(
    condition_id = NA,
    rep_id = NA,
    method = method,
    dx = dx,
    dy = dy,
    cell_size,
    RMSE = sqrt(mean(error ^ 2, na.rm = TRUE)),
    NRMSE_mean = sqrt(mean(error ^ 2, na.rm = TRUE)) / mean(y_true),
    NRMSE_range = sqrt(mean(error ^ 2, na.rm = TRUE)) / (max(y_true) - min(y_true)),
    NRMSE_sum = sqrt(mean(error ^ 2, na.rm = TRUE)) / sum(y_true),
    MAE = mean(abs(error), na.rm = TRUE),
    Mean_Error = mean(error, na.rm = TRUE),
    Median_Error = median(error, na.rm = TRUE),
    Min_Error = min(error),
    Max_Error = max(error),
    Skewness_Error = skewness(error, na.rm = T),
    Kurtosis_Error = kurtosis(error, na.rm = T),
    Variance_Error = var(error, na.rm = TRUE),
    Q1_Error = quantile(error, 0.25, na.rm = TRUE),
    Q3_Error = quantile(error, 0.75, na.rm = TRUE),
    R2 = cor(y_hat, y_true, use = "complete.obs") ^ 2,
    MAPE = mean(abs(error / y_true), na.rm = TRUE) * 100
  )
}
