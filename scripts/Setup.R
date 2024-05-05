# pipeline variables
start.time <- Sys.time()

# Set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path))

# since moving script from local to github - I want to adjust work dir to be main github dir - therefore 
setwd("..")

working.dir <- getwd()

# For reproducibility
set.seed(42)

# Explicitly load key packages
library("Seurat")
library("SeuratDisk")
library("dplyr")
library("scCustomize")
library("usefulfunctions")
library("ggplot2")
library("usefulfunctions")

# Establish colour scheme
source("scripts/variables/Colour_scheme_variable.R", local = knitr::knit_global())


# ggplot2 
# some issues with ggplot2 versions
# Nebulosa / scCustomize wont work with v3.5.0 and therefore needs 3.4.4
# many things can work with 3.4.4 version but
# clusttree requres ggplot2 >= 3.5 

# therefore function to change the version of ggplot2 as required so code can function
ggplot2_switch_version <- function(desired_version) {
  
  library("devtools")
  
  # Get the current version of ggplot2
  current_version <- packageVersion("ggplot2")
  
  # Convert the current version to character to compare with desired version
  current_version <- as.character(current_version)
  
  # Check if the current version is different from the desired version
  if (current_version != desired_version) {
    print(paste("Current ggplot2 version is", current_version, ". Installing version", desired_version))
    # Install the desired version of ggplot2
    install_version("ggplot2", version = desired_version, repos = "http://cran.us.r-project.org")
    print(paste("ggplot2 version", desired_version, "has been installed."))
  } else {
    print(paste("ggplot2 version", desired_version, "is already installed."))
  }
}

########################################
# Parallelization of seurat functions
########################################

parallelization <- TRUE


if(parallelization){
  library(future)
  # check the current active plan
  plan()
  
  # detect cores
  ncores <- parallel::detectCores()
  ncores <- ncores-5
  # set plan
  plan("multisession", workers = ncores)
  plan()
  
  # set max limit for global variables
  # note this will drastically increase RAM usage, be aware of machine limits
  # set value is in bites therefore 1000*1024^2 = 1GB Therefore set to 10GBs 
  options(future.globals.maxSize = 10000 * 1024^2, 
          future.rng.onMisuse="ignore") 
}


###############################
# create output directories
###############################

# Common directories 
if(!dir.exists("saves")){dir.create("saves", recursive = T)}

# scRNAseq directories
if(!dir.exists("results/")){dir.create("results/", recursive = T)}


