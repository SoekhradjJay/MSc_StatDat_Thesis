# PROJECT DESCRIPTION

This repository contains the final code for my thesis: 
Evaluating the Performance of Simple Areal Weighting and Dasymetric Mapping for Population Redistribution under Global and Local Spatial Autocorrelation

The thesis was written in collaboration with Statistics Netherlands (Centraal Bureau voor de Statistiek).

# FOLDER STRUCTURE

After cloning the repository, the following structure can be found:

- `data` - preprocessing and results
  - `grids` - aggregation grids stored per translation and cell size in subfolders
  - `LU_poly` - land use polygons simulated in thesis
  - `results` - output of simulation stored per simulated scenario (called condition)
  - `study zone` - polygons of municipalities in The Netherlands, from: https://www.pdok.nl/atom-downloadservices/-/article/bestuurlijke-gebieden
- `plots_and_analysis` - .RMD files with code for images in methodology and results section
  - `To_git.R` - R.proj file which can be saved locally and used as reference root (recommended)
- `0.GetPackages.R` - R script Installer of all used packages for the simulation
- `1.Formulas&Functions.R` - R script for all formulas and functions used in the simulation. Alter parameters here.
- `2.Prep_StudyZone&LU&AggregationGrids.R` - Read (your own) study zone, land use polygons, and aggregation grids. Optionally make them.
- `3.Preprocessing.R` - Preprocessing steps (kriging + buffer zones)
- `4. RunSimulation.R` - Having run scripts 0-3, you can run the simulation.