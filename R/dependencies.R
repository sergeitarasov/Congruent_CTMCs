# library(rphenoscate)
# library(corHMM)
# library(phytools)
# library(plyr)
# library(dplyr)
# library(diversitree)
# library(hisse)

# List of packages for session
.packages <- c(
  'corHMM',
  'phytools',
  'plyr',
  'dplyr',
  'hisse',
  'diversitree'
)

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)

