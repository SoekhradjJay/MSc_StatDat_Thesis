################################################################################

# ------------------------- 0. Install/Load Packages  --------------------------

# This script is used to install and load packages necessary for simulation
# Last revision: 15 October 2025 (unless stated otherwise in GitHub)

################################################################################

# Get necessary packages
packages <- c("arrow", "data.table", "dplyr", "gstat", 
              "moments", "mvtnorm",
              "parallel", "raster", "sf", "sp", "terra", 
              "tidyverse", "future", "future.apply")

# function: install packages if missing 
install_if_missing <- function(pkgs) {
  installed <- rownames(installed.packages())
  to_install <- setdiff(pkgs, installed)
  if (length(to_install) > 0) install.packages(to_install)
}

# Execute installation (if missing)
install_if_missing(packages)

# Execute loading libraries
lapply(packages, library, character.only = TRUE)